"""Read the pixels of a PNG screenshot (standard library only), and classify them by colour.

decode_png is copied from the lab's WebGL checker (cf-lab: .claude/skills/webgl-threejs-graphics/scripts/check-webgl.py), which proved it on real
headless-Chromium screenshots. 8-bit RGB or RGBA, not interlaced: what Chromium writes."""
import struct
import zlib
from collections import Counter


def decode_png(data):
    """Return (width, height, [(r, g, b), ...]). Raises ValueError for anything this reader does not support."""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("not a PNG")
    pos, ihdr, idat = 8, None, b""
    while pos + 8 <= len(data):
        n, kind = struct.unpack(">I4s", data[pos:pos + 8])
        body = data[pos + 8:pos + 8 + n]
        pos += 12 + n
        if kind == b"IHDR":
            ihdr = struct.unpack(">IIBBBBB", body)
        elif kind == b"IDAT":
            idat += body
        elif kind == b"IEND":
            break
    if ihdr is None:
        raise ValueError("no IHDR")
    width, height, depth, ctype, _comp, _filt, interlace = ihdr
    if depth != 8 or ctype not in (2, 6) or interlace != 0:
        raise ValueError(f"unsupported PNG (depth {depth}, colour type {ctype}, interlace {interlace})")
    bpp = 3 if ctype == 2 else 4
    raw = zlib.decompress(idat)
    stride = width * bpp
    if len(raw) != height * (stride + 1):
        raise ValueError("PNG data has the wrong length")
    prev = bytearray(stride)
    pixels = []
    for y in range(height):
        f = raw[y * (stride + 1)]
        row = bytearray(raw[y * (stride + 1) + 1:(y + 1) * (stride + 1)])
        for i in range(stride):
            a = row[i - bpp] if i >= bpp else 0
            b = prev[i]
            c = prev[i - bpp] if i >= bpp else 0
            if f == 1:
                row[i] = (row[i] + a) & 255
            elif f == 2:
                row[i] = (row[i] + b) & 255
            elif f == 3:
                row[i] = (row[i] + (a + b) // 2) & 255
            elif f == 4:
                p = a + b - c
                pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                row[i] = (row[i] + (a if pa <= pb and pa <= pc else (b if pb <= pc else c))) & 255
            elif f != 0:
                raise ValueError(f"bad PNG filter {f}")
        pixels += [(row[i], row[i + 1], row[i + 2]) for i in range(0, stride, bpp)]
        prev = row
    return width, height, pixels


def pixels(png):
    """(width, height, list of (r, g, b))."""
    return decode_png(png)


def nonblank(png, background_max=24):
    """Fraction of pixels that are not near-black."""
    w, h, px = decode_png(png)
    return sum(1 for r, g, b in px if max(r, g, b) > background_max) / (w * h)


def rose_points(png):
    """(x, y) of every saturated rose-coloured pixel (the residue 508 marker), as lists, for measuring the marker's shape."""
    import colorsys
    w, h, px = decode_png(png)
    pts = []
    for i, (r, g, b) in enumerate(px):
        hh, ss, vv = colorsys.rgb_to_hsv(r / 255, g / 255, b / 255)
        if vv >= 0.8 and ss >= 0.35 and (hh * 360 >= 335 or hh * 360 < 10):          # bright: the faint roses behind the transparent canvas are not the marker
            pts.append((i % w, i // w))
    return pts


def white_count(png, floor=235):
    """How many pixels are near-white: the walk marker is white, the structure is saturated colour."""
    w, h, px = decode_png(png)
    return sum(1 for r, g, b in px if min(r, g, b) >= floor)


def hue_counts(png):
    """How many saturated pixels fall in each colour family the viewer uses. By hue, because the viewer fades distant lines (brightness changes, hue does not)."""
    import colorsys
    w, h, px = decode_png(png)
    c = Counter()
    for r, g, b in px:
        hh, ss, vv = colorsys.rgb_to_hsv(r / 255, g / 255, b / 255)
        if vv < 0.3 or ss < 0.45:
            continue
        deg = hh * 360
        if 20 <= deg < 50:
            c["amber"] += 1
        elif deg >= 335 or deg < 10:
            c["rose"] += 1
        elif 165 <= deg < 197:
            c["cyan"] += 1
        elif 205 <= deg < 240:
            c["blue"] += 1
        elif 265 <= deg < 315:
            c["violet"] += 1
    return c
