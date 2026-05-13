"""Tests for genome_plotter.functions.utils."""

from __future__ import annotations

import unittest

from genome_plotter.functions.utils import Scaler


class TestScaler(unittest.TestCase):
    """Tests for Scaler.scale_genomic_coordinate."""

    def test_coordinate_at_gene_start_is_zero(self) -> None:
        """A coordinate equal to gene_start maps to pixel position 0."""
        s = Scaler(gene_start=1000, chunk_size=500, pixel_size=9)
        self.assertEqual(s.scale_genomic_coordinate(1000), 0)

    def test_one_full_chunk_equals_one_pixel_size(self) -> None:
        """One chunk ahead of gene_start maps to exactly one pixel_size."""
        # (1500 - 1000) / 500 * 9 = 9
        s = Scaler(gene_start=1000, chunk_size=500, pixel_size=9)
        self.assertEqual(s.scale_genomic_coordinate(1500), 9)

    def test_half_chunk_is_half_pixel(self) -> None:
        """Half a chunk ahead maps to half a pixel_size."""
        # 500 / 1000 * 10 = 5
        s = Scaler(gene_start=0, chunk_size=1000, pixel_size=10)
        self.assertEqual(s.scale_genomic_coordinate(500), 5)

    def test_coordinate_before_start_is_negative(self) -> None:
        """A coordinate before gene_start produces a negative pixel value."""
        # (500 - 1000) / 500 * 9 = -9
        s = Scaler(gene_start=1000, chunk_size=500, pixel_size=9)
        self.assertEqual(s.scale_genomic_coordinate(500), -9)

    def test_fractional_result_is_rounded(self) -> None:
        """Non-integer pixel values are rounded to the nearest integer."""
        # 1 / 3 * 2 = 0.666… → rounds to 1
        s = Scaler(gene_start=0, chunk_size=3, pixel_size=2)
        self.assertEqual(s.scale_genomic_coordinate(1), 1)

    def test_zero_gene_start(self) -> None:
        """Scaler works correctly when gene_start is 0."""
        s = Scaler(gene_start=0, chunk_size=450, pixel_size=18)
        self.assertEqual(s.scale_genomic_coordinate(450 * 200), 18 * 200)

    def test_result_is_int(self) -> None:
        """Return value is always a Python int."""
        s = Scaler(gene_start=0, chunk_size=100, pixel_size=5)
        self.assertIsInstance(s.scale_genomic_coordinate(50), int)


if __name__ == "__main__":
    unittest.main()
