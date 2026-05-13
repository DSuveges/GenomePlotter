"""Tests for genome_plotter.input_parsers.fetch_gwas_catalog."""

from __future__ import annotations

import unittest

import pandas as pd

from genome_plotter.input_parsers.fetch_gwas_catalog import FetchGwas


class _FakeGwas(FetchGwas):
    """Minimal FetchGwas subclass for testing process_gwas_data without FTP."""

    def __init__(self, raw: pd.DataFrame) -> None:
        """Initialise with raw data only; skip all FTP setup."""
        self.raw_data = raw


class TestFetchGwasProcessData(unittest.TestCase):
    """Tests for FetchGwas.process_gwas_data.

    process_gwas_data reads from self.raw_data and writes to self.gwas_df.
    _FakeGwas subclasses FetchGwas and overrides __init__ to bypass FTP.
    """

    def _run(self, rows: list[dict[str, object]]) -> pd.DataFrame:
        """Build a raw DataFrame, call process_gwas_data, return gwas_df."""
        obj = _FakeGwas(pd.DataFrame(rows).copy())
        obj.process_gwas_data()
        return obj.gwas_df

    def _valid(self, **overrides: object) -> dict[str, object]:
        """Return a valid raw GWAS row, optionally overriding fields."""
        base: dict[str, object] = {"CHR_ID": "1", "CHR_POS": 1000, "SNPS": "rs12345"}
        base.update(overrides)
        return base

    def test_valid_row_passes_through(self) -> None:
        """A fully valid row survives all filters."""
        result = self._run([self._valid()])
        self.assertEqual(len(result), 1)

    def test_output_columns(self) -> None:
        """Output DataFrame has exactly the expected four columns."""
        result = self._run([self._valid()])
        self.assertListEqual(list(result.columns), ["#chr", "start", "end", "rsID"])

    def test_end_equals_start_position(self) -> None:
        """start and end are both set to CHR_POS (single-base position)."""
        result = self._run([self._valid(CHR_POS=5000)])
        self.assertEqual(result.iloc[0]["start"], 5000)
        self.assertEqual(result.iloc[0]["end"], 5000)

    def test_null_chr_pos_is_filtered(self) -> None:
        """Rows with a null CHR_POS are removed."""
        result = self._run([self._valid(), self._valid(CHR_POS=None)])
        self.assertEqual(len(result), 1)

    def test_chr_id_with_lowercase_x_is_filtered(self) -> None:
        """CHR_ID containing lowercase 'x' (malformed entry) is removed."""
        result = self._run([self._valid(), self._valid(CHR_ID="x")])
        self.assertEqual(len(result), 1)

    def test_chr_id_with_semicolon_is_filtered(self) -> None:
        """CHR_ID containing ';' (multi-locus entry) is removed."""
        result = self._run([self._valid(), self._valid(CHR_ID="1;2")])
        self.assertEqual(len(result), 1)

    def test_snps_without_rs_prefix_is_filtered(self) -> None:
        """SNPS not containing 'rs' are removed."""
        result = self._run([self._valid(), self._valid(SNPS="chr1:12345")])
        self.assertEqual(len(result), 1)

    def test_snps_rs_match_is_case_insensitive(self) -> None:
        """'rs' match in SNPS is case-insensitive."""
        result = self._run([self._valid(SNPS="RS9876")])
        self.assertEqual(len(result), 1)

    def test_duplicate_rows_are_deduplicated(self) -> None:
        """Identical rows are collapsed to a single entry."""
        result = self._run([self._valid(), self._valid()])
        self.assertEqual(len(result), 1)

    def test_multiple_chromosomes_preserved(self) -> None:
        """Rows on different chromosomes are all retained."""
        result = self._run([self._valid(CHR_ID="1"), self._valid(CHR_ID="2")])
        self.assertEqual(len(result), 2)
        self.assertIn("1", result["#chr"].values)
        self.assertIn("2", result["#chr"].values)


if __name__ == "__main__":
    unittest.main()
