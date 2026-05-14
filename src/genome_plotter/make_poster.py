"""Compose per-chromosome PNGs into a two-row genome poster."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from loguru import logger
from PIL import Image, ImageDraw, ImageFont

from genome_plotter import LOG_FORMAT

# ---------------------------------------------------------------------------
# Chromosome layout
# ---------------------------------------------------------------------------

# Each column is a list of chromosome names.  Two names = stacked vertically.
_ROW1: list[list[str]] = [[str(i)] for i in range(1, 13)]       # chr1–chr12
_ROW2: list[list[str]] = [
    ["22"],                                                        # under chr1
    *[[str(i)] for i in range(13, 22)],                           # chr13–chr21
    ["X", "Y"],                                                    # stacked last column
]

# Gap constants are expressed in units of the *original* chromosome image width
# so that they scale proportionally with TARGET_CHR_WIDTH.
_H_GAP_FRAC = 0.027       # horizontal gap between columns  (~60 / 2241)
_STACK_GAP_FRAC = 0.013   # gap between stacked images      (~30 / 2241)
_V_GAP_FRAC = 0.062       # vertical gap between the two rows (~140 / 2241)
_MARGIN_FRAC = 0.036      # outer margin                     (~80 / 2241)

_LABEL_MARGIN_PX = 14     # fixed px below each image before its label
_FONT_SIZE_PT = 28        # fixed pt — not scaled so labels stay legible

_BACKGROUND = (255, 255, 255)
_LABEL_COLOR = (40, 40, 40)


# ---------------------------------------------------------------------------
# Font discovery
# ---------------------------------------------------------------------------

# Common TrueType font paths, tried in order across macOS / Linux / Windows.
_FONT_CANDIDATES = [
    # macOS
    "/System/Library/Fonts/SFNSMono.ttf",
    "/System/Library/Fonts/Geneva.ttf",
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    # Linux – Debian / Ubuntu
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/truetype/freefont/FreeSans.ttf",
    # Linux – Arch
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
    # Linux – Fedora / RHEL
    "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
    # Windows
    "C:/Windows/Fonts/arial.ttf",
    "C:/Windows/Fonts/calibri.ttf",
]


def find_font(
    size: int, path: str | None = None
) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Return a PIL font at *size* pt.

    If *path* is given it is tried first and an error is raised if it fails.
    Otherwise, a platform-appropriate system font is located automatically,
    falling back to PIL's built-in bitmap font if nothing else works.

    Args:
        size: Desired font size in points.
        path: Optional explicit path to a TrueType/OpenType font file.

    Returns:
        A PIL font object.

    Raises:
        ValueError: If *path* is given but the file cannot be loaded.
    """
    if path:
        try:
            return ImageFont.truetype(path, size)
        except (OSError, IOError) as exc:
            raise ValueError(f"Cannot load font from {path!r}: {exc}") from exc

    for candidate in _FONT_CANDIDATES:
        try:
            return ImageFont.truetype(candidate, size)
        except (OSError, IOError):
            continue

    logger.warning(
        "No TrueType font found on this system; falling back to PIL built-in bitmap font. "
        "Labels may look pixelated. Pass --font <path> to specify a font."
    )
    return ImageFont.load_default()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load(data_dir: Path, name: str) -> Image.Image | None:
    """Load a chromosome PNG, returning None if the file does not exist."""
    path = data_dir / f"chr{name}.png"
    if not path.exists():
        logger.warning(f"{path} not found – skipping.")
        return None
    return Image.open(path).convert("RGBA")


def _scale(img: Image.Image, factor: float) -> Image.Image:
    """Return a copy of *img* resized by *factor* using Lanczos resampling."""
    return img.resize(
        (max(1, round(img.width * factor)), max(1, round(img.height * factor))),
        Image.LANCZOS,
    )


def _text_width(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.FreeTypeFont | ImageFont.ImageFont) -> int:
    """Return the rendered pixel width of *text* with *font*."""
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0]


def _col_height(imgs: list[Image.Image], stack_gap: int) -> int:
    """Return total pixel height of a column including inter-image gaps."""
    return sum(img.height for img in imgs) + stack_gap * (len(imgs) - 1)


def _row_dims(
    layout: list[list[str]],
    loaded: dict[str, Image.Image],
    h_gap: int,
    stack_gap: int,
) -> tuple[int, int]:
    """Return (total_width, max_column_height) for a row."""
    widths, heights = [], []
    for col in layout:
        imgs = [loaded[n] for n in col if n in loaded]
        if imgs:
            widths.append(imgs[0].width)
            heights.append(_col_height(imgs, stack_gap))
    total_w = sum(widths) + h_gap * (len(widths) - 1)
    return total_w, max(heights, default=0)


def _draw_row(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    layout: list[list[str]],
    loaded: dict[str, Image.Image],
    row_top: int,
    row_max_h: int,
    align: str,
    h_gap: int,
    stack_gap: int,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    canvas_w: int,
) -> None:
    """Paste one row of chromosome images onto *canvas* and draw labels."""
    present = [
        ([n for n in col if n in loaded], [loaded[n] for n in col if n in loaded])
        for col in layout
        if any(n in loaded for n in col)
    ]

    total_w = sum(imgs[0].width for _, imgs in present) + h_gap * (len(present) - 1)
    x = (canvas_w - total_w) // 2

    for names, imgs in present:
        col_h = _col_height(imgs, stack_gap)
        col_w = imgs[0].width
        y = row_top if align == "top" else row_top + row_max_h - col_h

        for name, img in zip(names, imgs):
            canvas.paste(img, (x, y), img)
            label = f"chr{name}"
            lx = x + (col_w - _text_width(draw, label, font)) // 2
            draw.text((lx, y + img.height + _LABEL_MARGIN_PX), label, fill=_LABEL_COLOR, font=font)
            y += img.height + stack_gap

        x += col_w + h_gap


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_poster(
    data_dir: Path,
    target_chr_width: int = 380,
    font_path: str | None = None,
) -> Image.Image:
    """Compose a two-row genome poster from per-chromosome PNGs.

    Args:
        data_dir: Directory containing ``chr<N>.png`` files.
        target_chr_width: Desired pixel width of each chromosome column.
        font_path: Optional path to a TrueType font; auto-detected if omitted.

    Returns:
        A PIL ``Image`` containing the composed poster.

    Raises:
        FileNotFoundError: If no chromosome PNGs are found in *data_dir*.
    """
    all_names = {n for col in _ROW1 + _ROW2 for n in col}
    raw = {n: _load(data_dir, n) for n in all_names}
    sample = next((img for img in raw.values() if img is not None), None)
    if sample is None:
        raise FileNotFoundError(f"No chromosome PNGs found in {data_dir}")

    scale = target_chr_width / sample.width
    logger.info(f"Scale factor: {scale:.4f}  (chr width → {target_chr_width}px)")

    loaded = {n: _scale(img, scale) for n, img in raw.items() if img is not None}

    chr_w = next(iter(loaded.values())).width
    h_gap = max(1, round(_H_GAP_FRAC * chr_w))
    stack_gap = max(1, round(_STACK_GAP_FRAC * chr_w))
    v_gap = max(1, round(_V_GAP_FRAC * chr_w))
    margin = max(1, round(_MARGIN_FRAC * chr_w))

    font = find_font(_FONT_SIZE_PT, font_path)

    r1_w, r1_max_h = _row_dims(_ROW1, loaded, h_gap, stack_gap)
    r2_w, r2_max_h = _row_dims(_ROW2, loaded, h_gap, stack_gap)

    max_stack = max(len(col) for col in _ROW1 + _ROW2)
    label_area = (_FONT_SIZE_PT + _LABEL_MARGIN_PX) * max_stack

    poster_w = max(r1_w, r2_w) + 2 * margin
    poster_h = margin + r1_max_h + label_area + v_gap + r2_max_h + label_area + margin

    canvas = Image.new("RGB", (poster_w, poster_h), _BACKGROUND)
    draw = ImageDraw.Draw(canvas)

    _draw_row(canvas, draw, _ROW1, loaded, margin, r1_max_h, "top", h_gap, stack_gap, font, poster_w)

    row2_top = margin + r1_max_h + label_area + v_gap
    _draw_row(canvas, draw, _ROW2, loaded, row2_top, r2_max_h, "bottom", h_gap, stack_gap, font, poster_w)

    return canvas


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments for make-poster.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Compose per-chromosome PNGs into a two-row genome poster."
    )
    parser.add_argument(
        "-f", "--folder",
        help="Directory containing chr<N>.png files (same folder used by plot-chromosome).",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file path (default: <folder>/genome_poster.png).",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-w", "--width",
        help="Target pixel width of each chromosome column (default: 380).",
        type=int,
        default=380,
    )
    parser.add_argument(
        "--font",
        help="Path to a TrueType/OpenType font file for labels. "
             "Auto-detected from common system locations if omitted.",
        type=str,
        default=None,
    )
    return parser.parse_args()


def main() -> None:
    """Entry point for the make-poster CLI command."""
    logger.add("genome_plotter.log", level="DEBUG", format=LOG_FORMAT)

    args = parse_arguments()
    data_dir = Path(args.folder)
    output = Path(args.output) if args.output else data_dir / "genome_poster.png"

    if not data_dir.is_dir():
        logger.error(f"Folder not found: {data_dir}")
        sys.exit(1)

    logger.info(f"Building poster from PNGs in {data_dir}")
    poster = build_poster(data_dir, target_chr_width=args.width, font_path=args.font)
    poster.save(output, dpi=(150, 150))
    logger.info(f"Saved {output}  ({poster.width}×{poster.height} px)")


if __name__ == "__main__":
    main()
