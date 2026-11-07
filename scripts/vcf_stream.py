# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 netmaster296
"""Sequential BGZF reader for gnomAD sites VCFs.

Bio.bgzf is correct and is what scripts/vcf_query.py uses for one-region
seeks. A chromosome annotation pass must not seek once per variant. This
reader inflates each BGZF block with zlib and walks forward. The TBI linear
index is used only to jump a gap that contains no requested position.
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path

import vcf_query

_GZIP_MAGIC = b"\x1f\x8b"
_WINDOW = 1 << 14


class BgzfStream:
    """Forward BGZF reader. tell() is the file offset of the next block."""

    def __init__(self, path):
        self.path = Path(path)
        self._fh = open(self.path, "rb", buffering=1024 * 1024)
        self._buf = bytearray()
        self._pos = 0
        self._base = 0

    def close(self) -> None:
        self._fh.close()

    def tell(self) -> int:
        return self._base + self._pos

    def seek_file(self, offset: int) -> None:
        self._fh.seek(offset)
        self._buf.clear()
        self._pos = 0
        self._base = offset

    def _need(self, count: int) -> bool:
        while self._pos + count > len(self._buf):
            chunk = self._fh.read(4 << 20)
            if not chunk:
                return False
            if self._pos > (1 << 20):
                del self._buf[: self._pos]
                self._base += self._pos
                self._pos = 0
            self._buf.extend(chunk)
        return True

    def read_block(self) -> bytes | None:
        if not self._need(18):
            return None
        if self._buf[self._pos : self._pos + 2] != _GZIP_MAGIC:
            raise ValueError(f"not a BGZF block at {self.tell()} in {self.path}")
        flag = self._buf[self._pos + 3]
        if flag & 4 == 0:
            raise ValueError(f"BGZF block has no extra field at {self.tell()}")
        extra_len = struct.unpack_from("<H", self._buf, self._pos + 10)[0]
        header_len = 12 + extra_len
        if not self._need(header_len):
            raise ValueError(f"truncated BGZF header at {self.tell()}")
        bsize = None
        index = self._pos + 12
        end = self._pos + header_len
        while index + 4 <= end:
            si1, si2, slen = struct.unpack_from("<BBH", self._buf, index)
            if si1 == 66 and si2 == 67 and slen >= 2:
                bsize = struct.unpack_from("<H", self._buf, index + 4)[0]
            index += 4 + slen
        if bsize is None:
            raise ValueError(f"BGZF block has no BC field at {self.tell()}")
        total = bsize + 1
        if not self._need(total):
            raise ValueError(f"truncated BGZF block at {self.tell()}")
        start = self._pos + header_len
        compressed = bytes(self._buf[start : self._pos + total - 8])
        try:
            data = zlib.decompress(compressed, wbits=-15)
        except zlib.error as exc:
            raise ValueError(f"BGZF inflate failed at {self.tell()}") from exc
        self._pos += total
        return data


def iter_lines(path, start_virtual: int | None = None):
    """Yield uncompressed lines without the newline. A line may cross blocks."""
    stream = BgzfStream(path)
    try:
        pending = b""
        if start_virtual is not None:
            stream.seek_file(start_virtual >> 16)
            first = stream.read_block()
            if first is None:
                return
            offset = start_virtual & 0xFFFF
            if offset > len(first):
                raise ValueError(f"virtual offset {start_virtual} is past the BGZF block")
            pending = first[offset:]
        while True:
            block = stream.read_block()
            if block is None:
                if pending:
                    yield pending
                return
            if pending:
                block = pending + block
                pending = b""
            parts = block.split(b"\n")
            pending = parts.pop()
            for part in parts:
                yield part
    finally:
        stream.close()


def contig_lengths(path) -> dict[str, int]:
    """Contig lengths from ##contig header lines. Stops at the first data line."""
    lengths = {}
    for line in iter_lines(path):
        if line.startswith(b"##contig="):
            text = line.decode("ascii")
            ident = _header_attr(text, "ID")
            length = _header_attr(text, "length")
            if ident and length and length.isdigit():
                lengths[ident] = int(length)
            continue
        if not line.startswith(b"#"):
            break
    return lengths


def info_header(path) -> dict[str, dict[str, str]]:
    """INFO field id -> Number, Type, and Description. Stops at the first data line."""
    fields = {}
    for line in iter_lines(path):
        if line.startswith(b"##INFO=<"):
            text = line.decode("ascii", errors="replace")
            ident = _header_attr(text, "ID")
            if not ident:
                continue
            fields[ident] = {
                "Number": _header_attr(text, "Number") or "",
                "Type": _header_attr(text, "Type") or "",
                "Description": _header_attr(text, "Description") or "",
            }
            continue
        if not line.startswith(b"#"):
            break
    return fields


def _header_attr(text: str, name: str) -> str | None:
    token = name + "="
    start = text.find(token)
    if start < 0:
        return None
    start += len(token)
    if start < len(text) and text[start] == '"':
        end = text.find('"', start + 1)
        if end < 0:
            return text[start + 1 :]
        return text[start + 1 : end]
    end = start
    while end < len(text) and text[end] not in {",", ">"}:
        end += 1
    return text[start:end]


def vep_format(description: str) -> list[str]:
    """Subfield names from a VEP Description that contains 'Format: A|B|C'."""
    marker = "Format:"
    index = description.find(marker)
    if index < 0:
        return []
    text = description[index + len(marker) :].strip()
    text = text.rstrip(".").strip()
    if "|" not in text:
        return []
    return [part.strip() for part in text.split("|")]


def parse_info(text: str) -> dict[str, str]:
    """Exact INFO keys. AF_nfe and AF_nfe_XX stay distinct."""
    out = {}
    for field in text.split(";"):
        if not field or field == ".":
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            out[key] = value
        else:
            out[field] = "1"
    return out


def allele_value(raw: str | None, index: int, number: str) -> str | None:
    """Pick one allele from a Number=A field. A single value is not indexed."""
    if raw is None or raw == "" or raw == ".":
        return None
    if number == "A" and "," in raw:
        parts = raw.split(",")
        if index < 0 or index >= len(parts):
            return None
        value = parts[index]
    else:
        value = raw
    if value == "" or value == ".":
        return None
    return value


class _LineScan:
    """One forward pass. A large empty gap seeks through the linear index."""

    def __init__(self, path, chrom: str):
        self.path = path
        self.chrom = chrom
        self.encoded = (chrom + "\t").encode("ascii")
        self.stream = BgzfStream(path)
        self.pending = b""
        self.index = vcf_query.TabixIndex(Path(str(path) + ".tbi"))
        self._started = False
        self._queue: list[bytes] = []
        self._qi = 0

    def close(self) -> None:
        self.stream.close()

    def _queue_block(self, block: bytes) -> None:
        if self.pending:
            block = self.pending + block
            self.pending = b""
        parts = block.split(b"\n")
        self.pending = parts.pop()
        self._queue = parts
        self._qi = 0

    def seek_pos(self, pos: int) -> None:
        virtual = self.index.virtual_offset(self.chrom, pos)
        if virtual is None:
            return
        block_off = virtual >> 16
        if self._started and block_off <= self.stream.tell():
            return
        self.stream.seek_file(block_off)
        self.pending = b""
        self._queue = []
        self._qi = 0
        first = self.stream.read_block()
        self._started = True
        if first is None:
            return
        offset = virtual & 0xFFFF
        chunk = first[offset:]
        parts = chunk.split(b"\n")
        # A virtual offset is a line start. Drop a fragment if it is not.
        if offset and parts and not (
            parts[0].startswith(self.encoded) or parts[0].startswith(b"#")
        ):
            parts = parts[1:]
        if not parts:
            self.pending = b""
            return
        self.pending = parts.pop()
        self._queue = parts
        self._qi = 0

    def next_line(self) -> bytes | None:
        while True:
            if self._qi < len(self._queue):
                line = self._queue[self._qi]
                self._qi += 1
                if line.startswith(self.encoded):
                    return line
                continue
            block = self.stream.read_block()
            if block is None:
                if self.pending.startswith(self.encoded):
                    line = self.pending
                    self.pending = b""
                    return line
                return None
            self._queue_block(block)


def stream_positions(path, chrom: str, positions: list[int], log_every: int = 2_000_000):
    """Yield (pos, line) for data lines whose POS is in the sorted position list.

    One forward pass. A gap of at least 256 kb seeks once to that next window.
    It does not seek once per variant and it does not seek backward.
    """
    if not positions:
        return
    scan = _LineScan(path, chrom)
    try:
        scan.seek_pos(positions[0])
        cursor = 0
        last_seek_pos = positions[0]
        seen = 0
        while True:
            line = scan.next_line()
            if line is None:
                return
            seen += 1
            tab = line.find(b"\t", len(scan.encoded))
            if tab < 0:
                continue
            pos = int(line[len(scan.encoded) : tab])
            if log_every and seen % log_every == 0:
                print(f"  {chrom} scanned {seen} lines, at {pos}", flush=True)
            while cursor < len(positions) and positions[cursor] < pos:
                cursor += 1
            if cursor >= len(positions):
                return
            if pos == positions[cursor]:
                yield pos, line
                continue
            nxt = positions[cursor]
            if nxt - pos >= 256_000 and nxt > last_seek_pos:
                scan.seek_pos(nxt)
                last_seek_pos = nxt
    finally:
        scan.close()
