"""Tests for genome_plotter.input_parsers.data_integrator."""

from __future__ import annotations

import unittest

import pandas as pd

from genome_plotter.input_parsers.data_integrator import DataIntegrator

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _genome_df(
    n: int = 10,
    chromosome: str = "1",
    chunk_size: int = 1000,
    gc_ratio: float | None = 0.5,
) -> pd.DataFrame:
    """Minimal genome DataFrame with n consecutive chunks."""
    starts = [i * chunk_size for i in range(n)]
    ends = [(i + 1) * chunk_size for i in range(n)]
    return pd.DataFrame(
        {
            "chr": chromosome,
            "start": starts,
            "end": ends,
            "GC_ratio": [gc_ratio] * n,
        }
    )


def _cytoband_df(
    chromosome: str = "1",
    centromere_start: int = 3000,
    centromere_end: int = 5000,
) -> pd.DataFrame:
    """Minimal cytoband DataFrame with one centromere spanning centromere_start–centromere_end."""
    mid = (centromere_start + centromere_end) // 2
    return pd.DataFrame(
        {
            "chr": [chromosome] * 4,
            "start": [0, 1000, centromere_start, mid],
            "end": [1000, centromere_start, mid, centromere_end],
            "name": ["p1", "p2", "p11", "q11"],
            "type": ["gneg", "gpos25", "acen", "acen"],
        }
    )


# ---------------------------------------------------------------------------
# DataIntegrator.add_xy_coordinates
# ---------------------------------------------------------------------------

class TestAddXYCoordinates(unittest.TestCase):
    """Tests for DataIntegrator.add_xy_coordinates."""

    def setUp(self) -> None:
        """Create a 6-chunk integrator for use across tests."""
        self.integrator = DataIntegrator(_genome_df(n=6))

    def test_first_chunk_is_at_origin(self) -> None:
        """The first chunk maps to x=0, y=0."""
        self.integrator.add_xy_coordinates(width=3)
        row = self.integrator.get_data().iloc[0]
        self.assertEqual(row["x"], 0)
        self.assertEqual(row["y"], 0)

    def test_x_wraps_at_width(self) -> None:
        """x resets to 0 at every multiple of width."""
        self.integrator.add_xy_coordinates(width=3)
        # index 3 → 3 % 3 = 0
        self.assertEqual(self.integrator.get_data().iloc[3]["x"], 0)

    def test_y_increments_at_row_boundary(self) -> None:
        """y increments by 1 each time x wraps around."""
        self.integrator.add_xy_coordinates(width=3)
        # index 3 → 3 // 3 = 1
        self.assertEqual(self.integrator.get_data().iloc[3]["y"], 1)

    def test_last_chunk_of_first_row(self) -> None:
        """Last chunk in the first row has x=width-1 and y=0."""
        self.integrator.add_xy_coordinates(width=3)
        row = self.integrator.get_data().iloc[2]
        self.assertEqual(row["x"], 2)
        self.assertEqual(row["y"], 0)

    def test_width_zero_puts_all_chunks_in_one_row(self) -> None:
        """width=0 treats the entire chromosome as one row (y=0 for all)."""
        self.integrator.add_xy_coordinates(width=0)
        self.assertTrue((self.integrator.get_data()["y"] == 0).all())

    def test_recompute_drops_old_coordinates(self) -> None:
        """Calling add_xy_coordinates again replaces previous x/y values."""
        self.integrator.add_xy_coordinates(width=3)
        self.integrator.add_xy_coordinates(width=2)
        # With width=2: index 2 → x=0, y=1
        row = self.integrator.get_data().iloc[2]
        self.assertEqual(row["x"], 0)
        self.assertEqual(row["y"], 1)

    def test_y_column_is_integer_type(self) -> None:
        """y column is stored as an integer dtype."""
        self.integrator.add_xy_coordinates(width=3)
        self.assertTrue(pd.api.types.is_integer_dtype(self.integrator.get_data()["y"]))


# ---------------------------------------------------------------------------
# DataIntegrator.add_centromere
# ---------------------------------------------------------------------------

class TestAddCentromere(unittest.TestCase):
    """Tests for DataIntegrator.add_centromere.

    Genome: 10 chunks × 1000 bp (0–10000).
    Centromere: 3000–5000 → overlapping chunks at index 3 (3000–4000) and 4 (4000–5000).
    Overlap condition: chunk.end > centromere_start AND chunk.start < centromere_end.
    """

    def setUp(self) -> None:
        """Create a 10-chunk integrator with a 3000–5000 centromere."""
        self.integrator = DataIntegrator(_genome_df(n=10))
        self.cytoband = _cytoband_df(centromere_start=3000, centromere_end=5000)

    def test_chunks_inside_centromere_are_marked(self) -> None:
        """Chunks overlapping the centromere region receive GENCODE='centromere'."""
        self.integrator.add_centromere(self.cytoband)
        data = self.integrator.get_data()
        self.assertEqual(data.iloc[3]["GENCODE"], "centromere")
        self.assertEqual(data.iloc[4]["GENCODE"], "centromere")

    def test_chunks_before_centromere_are_not_marked(self) -> None:
        """Chunks entirely before the centromere are not marked."""
        self.integrator.add_centromere(self.cytoband)
        data = self.integrator.get_data()
        for idx in [0, 1, 2]:
            self.assertNotEqual(data.iloc[idx]["GENCODE"], "centromere")

    def test_chunks_after_centromere_are_not_marked(self) -> None:
        """Chunks entirely after the centromere are not marked."""
        self.integrator.add_centromere(self.cytoband)
        data = self.integrator.get_data()
        for idx in [5, 6, 7, 8, 9]:
            self.assertNotEqual(data.iloc[idx]["GENCODE"], "centromere")

    def test_gencode_column_is_created_when_absent(self) -> None:
        """GENCODE column is created if it did not exist before."""
        self.assertNotIn("GENCODE", self.integrator.get_data().columns)
        self.integrator.add_centromere(self.cytoband)
        self.assertIn("GENCODE", self.integrator.get_data().columns)

    def test_centromere_overwrites_existing_gencode(self) -> None:
        """Centromere annotation overwrites any pre-existing GENCODE value."""
        self.integrator.__genome__["GENCODE"] = "intergenic"
        self.integrator.add_centromere(self.cytoband)
        data = self.integrator.get_data()
        self.assertEqual(data.iloc[3]["GENCODE"], "centromere")
        self.assertEqual(data.iloc[0]["GENCODE"], "intergenic")

    def test_only_matching_chromosome_is_used(self) -> None:
        """Centromere entries from other chromosomes are ignored."""
        extra = pd.concat([
            self.cytoband,
            pd.DataFrame({
                "chr": ["2", "2"],
                "start": [0, 500],
                "end": [500, 1000],
                "name": ["p11", "q11"],
                "type": ["acen", "acen"],
            }),
        ], ignore_index=True)
        self.integrator.add_centromere(extra)
        data = self.integrator.get_data()
        self.assertEqual(data.iloc[3]["GENCODE"], "centromere")
        self.assertNotEqual(data.iloc[0]["GENCODE"], "centromere")


# ---------------------------------------------------------------------------
# DataIntegrator.assign_hetero
# ---------------------------------------------------------------------------

class TestAssignHetero(unittest.TestCase):
    """Tests for DataIntegrator.assign_hetero."""

    def _integrator_with_gencode(
        self, n: int = 3, gc_ratios: list[float | None] | None = None
    ) -> DataIntegrator:
        """Return an integrator with a pre-populated GENCODE='intergenic' column."""
        df = _genome_df(n=n)
        if gc_ratios is not None:
            df["GC_ratio"] = gc_ratios
        integrator = DataIntegrator(df)
        integrator.__genome__["GENCODE"] = "intergenic"
        return integrator

    def test_null_gc_ratio_becomes_heterochromatin(self) -> None:
        """Chunks with null GC_ratio are labelled heterochromatin."""
        integrator = self._integrator_with_gencode(gc_ratios=[0.5, None, 0.5])
        integrator.assign_hetero()
        self.assertEqual(integrator.get_data().iloc[1]["GENCODE"], "heterochromatin")

    def test_non_null_gc_ratio_is_unchanged(self) -> None:
        """Chunks with a valid GC_ratio keep their existing GENCODE label."""
        integrator = self._integrator_with_gencode(gc_ratios=[0.5, None, 0.5])
        integrator.assign_hetero()
        data = integrator.get_data()
        self.assertEqual(data.iloc[0]["GENCODE"], "intergenic")
        self.assertEqual(data.iloc[2]["GENCODE"], "intergenic")

    def test_all_null_gc_all_become_heterochromatin(self) -> None:
        """When every chunk has null GC_ratio, all become heterochromatin."""
        integrator = self._integrator_with_gencode(n=4, gc_ratios=[None] * 4)
        integrator.assign_hetero()
        self.assertTrue((integrator.get_data()["GENCODE"] == "heterochromatin").all())

    def test_no_null_gc_nothing_changes(self) -> None:
        """When no chunk has null GC_ratio, no labels are changed."""
        integrator = self._integrator_with_gencode(n=4, gc_ratios=[0.4] * 4)
        integrator.assign_hetero()
        self.assertTrue((integrator.get_data()["GENCODE"] == "intergenic").all())


# ---------------------------------------------------------------------------
# DataIntegrator.add_dummy
# ---------------------------------------------------------------------------

class TestAddDummy(unittest.TestCase):
    """Tests for DataIntegrator.add_dummy."""

    def test_creates_gencode_column_when_absent(self) -> None:
        """GENCODE column is created if it did not exist."""
        integrator = DataIntegrator(_genome_df(n=3))
        self.assertNotIn("GENCODE", integrator.get_data().columns)
        integrator.add_dummy()
        self.assertIn("GENCODE", integrator.get_data().columns)

    def test_all_chunks_get_dummy_when_column_absent(self) -> None:
        """All chunks receive GENCODE='dummy' when no column existed."""
        integrator = DataIntegrator(_genome_df(n=5))
        integrator.add_dummy()
        self.assertTrue((integrator.get_data()["GENCODE"] == "dummy").all())

    def test_nan_values_filled_when_column_exists(self) -> None:
        """Null GENCODE values are replaced with 'dummy'."""
        integrator = DataIntegrator(_genome_df(n=3))
        integrator.__genome__["GENCODE"] = [None, "centromere", None]
        integrator.add_dummy()
        data = integrator.get_data()
        self.assertEqual(data.iloc[0]["GENCODE"], "dummy")
        self.assertEqual(data.iloc[2]["GENCODE"], "dummy")

    def test_existing_non_null_values_preserved(self) -> None:
        """Non-null GENCODE values are not overwritten by add_dummy."""
        integrator = DataIntegrator(_genome_df(n=3))
        integrator.__genome__["GENCODE"] = [None, "centromere", None]
        integrator.add_dummy()
        self.assertEqual(integrator.get_data().iloc[1]["GENCODE"], "centromere")

    def test_fully_populated_gencode_unchanged(self) -> None:
        """When every chunk already has a GENCODE label, nothing changes."""
        integrator = DataIntegrator(_genome_df(n=3))
        integrator.__genome__["GENCODE"] = ["exon", "centromere", "intergenic"]
        integrator.add_dummy()
        self.assertListEqual(
            list(integrator.get_data()["GENCODE"]),
            ["exon", "centromere", "intergenic"],
        )


if __name__ == "__main__":
    unittest.main()
