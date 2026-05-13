"""Tests for genome_plotter.functions.svg_handler."""

from __future__ import annotations

import unittest

from genome_plotter.functions.svg_handler import SvgHandler

# ---------------------------------------------------------------------------
# draw_rectangle
# ---------------------------------------------------------------------------

class TestDrawRectangle(unittest.TestCase):
    """Tests for SvgHandler.draw_rectangle."""

    def setUp(self) -> None:
        """Create a blank SvgHandler for each test."""
        self.svg = SvgHandler("", 500, 500)

    def test_rect_tag_present(self) -> None:
        """A <rect> element is appended to the SVG content."""
        self.svg.draw_rectangle(10, 20, 100, 50, "#000000", "#FF0000")
        self.assertIn("<rect", self.svg.get_svg())

    def test_coordinates_in_output(self) -> None:
        """x and y coordinates appear as SVG attributes."""
        self.svg.draw_rectangle(10, 20, 100, 50, "#000000", "#FF0000")
        content = self.svg.get_svg()
        self.assertIn('x="10"', content)
        self.assertIn('y="20"', content)

    def test_dimensions_in_output(self) -> None:
        """width and height appear as SVG attributes."""
        self.svg.draw_rectangle(10, 20, 100, 50, "#000000", "#FF0000")
        content = self.svg.get_svg()
        self.assertIn('width="100"', content)
        self.assertIn('height="50"', content)

    def test_stroke_color_in_output(self) -> None:
        """Stroke color is embedded in the style attribute."""
        self.svg.draw_rectangle(0, 0, 10, 10, "#ABCDEF", "#000000")
        self.assertIn("#ABCDEF", self.svg.get_svg())

    def test_fill_color_in_output(self) -> None:
        """Fill color is embedded in the style attribute."""
        self.svg.draw_rectangle(0, 0, 10, 10, "#000000", "#123456")
        self.assertIn("#123456", self.svg.get_svg())

    def test_multiple_rects_accumulate(self) -> None:
        """Each draw_rectangle call appends a new <rect> element."""
        self.svg.draw_rectangle(0, 0, 10, 10, "#000", "#FFF")
        self.svg.draw_rectangle(20, 20, 30, 30, "#000", "#FFF")
        self.assertEqual(self.svg.get_svg().count("<rect"), 2)


# ---------------------------------------------------------------------------
# draw_line
# ---------------------------------------------------------------------------

class TestDrawLine(unittest.TestCase):
    """Tests for SvgHandler.draw_line."""

    def setUp(self) -> None:
        """Create a blank SvgHandler for each test."""
        self.svg = SvgHandler("", 500, 500)

    def test_line_tag_present(self) -> None:
        """A <line> element is appended to the SVG content."""
        self.svg.draw_line(0, 0, 100, 100)
        self.assertIn("<line", self.svg.get_svg())

    def test_coordinates_in_output(self) -> None:
        """x1, y1, x2, y2 appear as SVG attributes."""
        self.svg.draw_line(10, 20, 30, 40)
        content = self.svg.get_svg()
        self.assertIn('x1="10"', content)
        self.assertIn('y1="20"', content)
        self.assertIn('x2="30"', content)
        self.assertIn('y2="40"', content)

    def test_default_stroke_color(self) -> None:
        """stroke defaults to #000000 when not specified."""
        self.svg.draw_line(0, 0, 10, 10)
        self.assertIn('stroke="#000000"', self.svg.get_svg())

    def test_custom_stroke_color(self) -> None:
        """A custom stroke color is written to the SVG attribute."""
        self.svg.draw_line(0, 0, 10, 10, stroke="#FF0000")
        self.assertIn('stroke="#FF0000"', self.svg.get_svg())

    def test_default_stroke_width(self) -> None:
        """stroke-width defaults to 3 when not specified."""
        self.svg.draw_line(0, 0, 10, 10)
        self.assertIn('stroke-width="3"', self.svg.get_svg())

    def test_custom_stroke_width(self) -> None:
        """A custom stroke_width is written to the SVG attribute."""
        self.svg.draw_line(0, 0, 10, 10, stroke_width=7)
        self.assertIn('stroke-width="7"', self.svg.get_svg())

    def test_kwarg_underscore_converted_to_hyphen(self) -> None:
        """Kwarg names have underscores replaced with hyphens in the SVG."""
        self.svg.draw_line(0, 0, 10, 10, stroke_dasharray="4")
        self.assertIn('stroke-dasharray="4"', self.svg.get_svg())

    def test_multiple_kwargs_all_converted(self) -> None:
        """Multiple kwargs are all appended with hyphenated names."""
        self.svg.draw_line(0, 0, 10, 10, stroke_dasharray="4", stroke_opacity="0.5")
        content = self.svg.get_svg()
        self.assertIn('stroke-dasharray="4"', content)
        self.assertIn('stroke-opacity="0.5"', content)

    def test_no_kwargs_leaves_no_stray_attributes(self) -> None:
        """When no kwargs are passed, no extra attribute fragments appear."""
        self.svg.draw_line(0, 0, 10, 10)
        self.assertNotIn("stroke-dasharray", self.svg.get_svg())


# ---------------------------------------------------------------------------
# add_text
# ---------------------------------------------------------------------------

class TestAddText(unittest.TestCase):
    """Tests for SvgHandler.add_text."""

    def setUp(self) -> None:
        """Create a blank SvgHandler for each test."""
        self.svg = SvgHandler("", 500, 500)

    def test_text_tag_present(self) -> None:
        """A <text> element is appended to the SVG content."""
        self.svg.add_text(50, 100, "Hello")
        self.assertIn("<text", self.svg.get_svg())

    def test_text_content_in_output(self) -> None:
        """The text string appears inside the <text> element."""
        self.svg.add_text(50, 100, "GenomePlotter")
        self.assertIn("GenomePlotter", self.svg.get_svg())

    def test_coordinates_in_output(self) -> None:
        """x and y coordinates appear as SVG attributes."""
        self.svg.add_text(50, 100, "x")
        content = self.svg.get_svg()
        self.assertIn('x="50"', content)
        self.assertIn('y="100"', content)

    def test_default_anchor_is_start(self) -> None:
        """text-anchor defaults to 'start' when not specified."""
        self.svg.add_text(0, 0, "x")
        self.assertIn('text-anchor="start"', self.svg.get_svg())

    def test_custom_anchor(self) -> None:
        """A custom anchor value is written to the text-anchor attribute."""
        self.svg.add_text(0, 0, "x", anchor="middle")
        self.assertIn('text-anchor="middle"', self.svg.get_svg())

    def test_font_size_in_output(self) -> None:
        """Font size is written as font-size in pixels."""
        self.svg.add_text(0, 0, "x", size=14)
        self.assertIn('font-size="14px"', self.svg.get_svg())

    def test_fill_color_in_output(self) -> None:
        """A custom fill color is written to the fill attribute."""
        self.svg.add_text(0, 0, "x", fill="#FF0000")
        self.assertIn('fill="#FF0000"', self.svg.get_svg())

    def test_default_fill_is_black(self) -> None:
        """fill defaults to #000000 when not specified."""
        self.svg.add_text(0, 0, "x")
        self.assertIn('fill="#000000"', self.svg.get_svg())


# ---------------------------------------------------------------------------
# append_svg
# ---------------------------------------------------------------------------

class TestAppendSvg(unittest.TestCase):
    """Tests for SvgHandler.append_svg."""

    def test_content_is_appended(self) -> None:
        """Appended string is joined to the existing SVG content."""
        svg = SvgHandler("initial", 100, 100)
        svg.append_svg(" appended")
        self.assertEqual(svg.get_svg(), "initial appended")

    def test_multiple_appends_accumulate(self) -> None:
        """Successive appends all accumulate in order."""
        svg = SvgHandler("", 100, 100)
        svg.append_svg("A")
        svg.append_svg("B")
        svg.append_svg("C")
        self.assertEqual(svg.get_svg(), "ABC")

    def test_dimensions_unchanged_by_append(self) -> None:
        """append_svg does not alter the canvas dimensions."""
        svg = SvgHandler("", 100, 200)
        svg.append_svg("<rect />")
        self.assertEqual(svg.get_width(), 100)
        self.assertEqual(svg.get_height(), 200)


# ---------------------------------------------------------------------------
# group
# ---------------------------------------------------------------------------

class TestGroup(unittest.TestCase):
    """Tests for SvgHandler.group."""

    def test_content_wrapped_in_g_transform(self) -> None:
        """Existing content is wrapped in a <g transform='translate(...)'> element."""
        svg = SvgHandler("content", 100, 200)
        svg.group(translate=(30, 10))
        self.assertIn('<g transform="translate(30 10)">', svg.get_svg())

    def test_original_content_preserved_inside_group(self) -> None:
        """The original SVG content appears inside the group element."""
        svg = SvgHandler("content", 100, 200)
        svg.group(translate=(30, 10))
        self.assertIn("content", svg.get_svg())

    def test_width_increased_by_x_translation(self) -> None:
        """Canvas width grows by the absolute x translation value."""
        svg = SvgHandler("", 100, 200)
        svg.group(translate=(50, 0))
        self.assertEqual(svg.get_width(), 150)

    def test_height_increased_by_y_translation(self) -> None:
        """Canvas height grows by the absolute y translation value."""
        svg = SvgHandler("", 100, 200)
        svg.group(translate=(0, 30))
        self.assertEqual(svg.get_height(), 230)

    def test_negative_translation_expands_canvas_via_abs(self) -> None:
        """Negative translation still expands the canvas (absolute value is used)."""
        svg = SvgHandler("", 100, 200)
        svg.group(translate=(-50, 0))
        self.assertEqual(svg.get_width(), 150)

    def test_zero_translation_leaves_dimensions_unchanged(self) -> None:
        """A zero translation does not alter the canvas dimensions."""
        svg = SvgHandler("", 100, 200)
        svg.group(translate=(0, 0))
        self.assertEqual(svg.get_width(), 100)
        self.assertEqual(svg.get_height(), 200)

    def test_default_translate_is_zero(self) -> None:
        """Calling group() without arguments applies a (0 0) translation."""
        svg = SvgHandler("content", 100, 200)
        svg.group()
        self.assertIn("translate(0 0)", svg.get_svg())


# ---------------------------------------------------------------------------
# merge_svg
# ---------------------------------------------------------------------------

class TestMergeSvg(unittest.TestCase):
    """Tests for SvgHandler.merge_svg."""

    def test_other_content_is_merged_into_self(self) -> None:
        """Content from both handlers appears in the merged result."""
        a = SvgHandler("AAA", 100, 100)
        b = SvgHandler("BBB", 50, 50)
        a.merge_svg(b)
        self.assertIn("AAA", a.get_svg())
        self.assertIn("BBB", a.get_svg())

    def test_width_takes_maximum(self) -> None:
        """Canvas width after merge is the maximum of the two widths."""
        a = SvgHandler("", 100, 100)
        b = SvgHandler("", 200, 50)
        a.merge_svg(b)
        self.assertEqual(a.get_width(), 200)

    def test_height_takes_maximum(self) -> None:
        """Canvas height after merge is the maximum of the two heights."""
        a = SvgHandler("", 100, 100)
        b = SvgHandler("", 50, 300)
        a.merge_svg(b)
        self.assertEqual(a.get_height(), 300)

    def test_self_wins_when_larger(self) -> None:
        """Self dimensions are preserved when larger than the other handler."""
        a = SvgHandler("", 300, 400)
        b = SvgHandler("", 100, 100)
        a.merge_svg(b)
        self.assertEqual(a.get_width(), 300)
        self.assertEqual(a.get_height(), 400)

    def test_other_is_not_mutated(self) -> None:
        """The merged-in handler is not altered by the operation."""
        a = SvgHandler("A", 100, 100)
        b = SvgHandler("B", 50, 50)
        a.merge_svg(b)
        self.assertEqual(b.get_svg(), "B")
        self.assertEqual(b.get_width(), 50)
        self.assertEqual(b.get_height(), 50)


# ---------------------------------------------------------------------------
# _build_svg
# ---------------------------------------------------------------------------

class TestBuildSvg(unittest.TestCase):
    """Tests for SvgHandler._build_svg."""

    def test_output_starts_with_svg_tag(self) -> None:
        """The returned string begins with an <svg ...> opening tag."""
        result = SvgHandler("", 100, 200)._build_svg()
        self.assertTrue(result.strip().startswith("<svg"))

    def test_output_ends_with_closing_svg_tag(self) -> None:
        """The returned string ends with a </svg> closing tag."""
        result = SvgHandler("", 100, 200)._build_svg()
        self.assertTrue(result.strip().endswith("</svg>"))

    def test_width_in_svg_tag(self) -> None:
        """Canvas width is written into the <svg> tag."""
        result = SvgHandler("", 123, 456)._build_svg()
        self.assertIn('width="123"', result)

    def test_height_in_svg_tag(self) -> None:
        """Canvas height is written into the <svg> tag."""
        result = SvgHandler("", 123, 456)._build_svg()
        self.assertIn('height="456"', result)

    def test_content_present_in_output(self) -> None:
        """SVG content passed at construction appears in the document."""
        result = SvgHandler("<rect />", 100, 100)._build_svg()
        self.assertIn("<rect />", result)

    def test_xmlns_in_svg_header(self) -> None:
        """The SVG namespace declaration is present in the header."""
        result = SvgHandler("", 100, 100)._build_svg()
        self.assertIn('xmlns="http://www.w3.org/2000/svg"', result)

    def test_no_background_rect_when_background_is_none(self) -> None:
        """No background <rect> is added when background=None."""
        result = SvgHandler("", 100, 100, background=None)._build_svg()
        self.assertNotIn('width="100%"', result)

    def test_background_fill_present_when_set(self) -> None:
        """The background color appears in a full-canvas <rect> fill."""
        result = SvgHandler("", 100, 100, background="#FFFFFF")._build_svg()
        self.assertIn('fill="#FFFFFF"', result)

    def test_background_rect_covers_full_canvas(self) -> None:
        """The background <rect> spans 100% width and height."""
        result = SvgHandler("", 100, 100, background="#FFFFFF")._build_svg()
        self.assertIn('width="100%"', result)
        self.assertIn('height="100%"', result)

    def test_content_appears_after_svg_opening_tag(self) -> None:
        """SVG content is positioned after the opening <svg> tag."""
        result = SvgHandler("MARKER", 100, 100)._build_svg()
        svg_open_pos = result.index("<svg")
        marker_pos = result.index("MARKER")
        self.assertGreater(marker_pos, svg_open_pos)


if __name__ == "__main__":
    unittest.main()
