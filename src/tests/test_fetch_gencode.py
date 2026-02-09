"""Tests for the FetchGencode data processing methods."""

from __future__ import annotations

import os
import unittest

import pandas as pd

from genome_plotter.input_parsers.fetch_gencode import FetchGencode

SAMPLE_GTF = os.path.join(
    os.path.dirname(__file__), "sample_data", "gencode_sample.gtf.gz"
)


class TestParseGtfAnnotations(unittest.TestCase):
    """Test cases for FetchGencode.parse_gtf_annotations."""

    def test_basic_parsing(self) -> None:
        """Annotation strings are split into separate columns."""
        df = pd.DataFrame(
            {
                "chr": ["1", "1"],
                "type": ["gene", "transcript"],
                "start": [100, 100],
                "end": [200, 200],
                "strand": ["+", "+"],
                "annotation": [
                    'gene_id "ENSG00000001.1"; gene_name "ABC"; gene_type "protein_coding"',
                    'gene_id "ENSG00000001.1"; gene_name "ABC"; transcript_id "ENST001"',
                ],
            }
        )
        result = FetchGencode.parse_gtf_annotations(df)

        self.assertNotIn("annotation", result.columns)
        self.assertIn("gene_id", result.columns)
        self.assertIn("gene_name", result.columns)
        self.assertEqual(result.iloc[0]["gene_id"], "ENSG00000001.1")
        self.assertEqual(result.iloc[0]["gene_name"], "ABC")
        self.assertEqual(result.iloc[0]["gene_type"], "protein_coding")

    def test_preserves_coordinate_columns(self) -> None:
        """Coordinate columns are preserved after parsing."""
        df = pd.DataFrame(
            {
                "chr": ["1"],
                "type": ["gene"],
                "start": [100],
                "end": [200],
                "strand": ["+"],
                "annotation": ['gene_id "ENSG001"; gene_name "X"'],
            }
        )
        result = FetchGencode.parse_gtf_annotations(df)
        for col in ["chr", "type", "start", "end", "strand"]:
            self.assertIn(col, result.columns)


class TestStripGeneIdVersion(unittest.TestCase):
    """Test cases for FetchGencode.strip_gene_id_version."""

    def test_version_stripped(self) -> None:
        """Version suffixes are removed from gene_id."""
        df = pd.DataFrame({"gene_id": ["ENSG00000001.5", "ENSG00000002.12"]})
        result = FetchGencode.strip_gene_id_version(df)
        self.assertEqual(result.gene_id.tolist(), ["ENSG00000001", "ENSG00000002"])

    def test_no_version(self) -> None:
        """gene_id without version suffix is unchanged."""
        df = pd.DataFrame({"gene_id": ["ENSG00000001"]})
        result = FetchGencode.strip_gene_id_version(df)
        self.assertEqual(result.gene_id.tolist(), ["ENSG00000001"])


class TestFilterProteinCodingGenes(unittest.TestCase):
    """Test cases for FetchGencode.filter_protein_coding_genes."""

    def _make_df(self, chr: str, gene_type: str, gene_name: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "chr": [chr],
                "gene_type": [gene_type],
                "gene_name": [gene_name],
                "type": ["gene"],
            }
        )

    def test_keeps_valid_row(self) -> None:
        """Standard protein-coding gene on conventional chromosome passes."""
        df = self._make_df("1", "protein_coding", "BRCA1")
        result = FetchGencode.filter_protein_coding_genes(df)
        self.assertEqual(len(result), 1)

    def test_filters_non_conventional_chr(self) -> None:
        """Entries on non-conventional chromosomes are removed."""
        df = self._make_df("KI270757.1", "protein_coding", "BRCA1")
        result = FetchGencode.filter_protein_coding_genes(df)
        self.assertEqual(len(result), 0)

    def test_filters_non_protein_coding(self) -> None:
        """Non-protein_coding gene_type entries are removed."""
        df = self._make_df("1", "lncRNA", "XIST")
        result = FetchGencode.filter_protein_coding_genes(df)
        self.assertEqual(len(result), 0)

    def test_filters_ensg_gene_name(self) -> None:
        """Gene names starting with ENSG are removed."""
        df = self._make_df("1", "protein_coding", "ENSG00000123456")
        result = FetchGencode.filter_protein_coding_genes(df)
        self.assertEqual(len(result), 0)

    def test_filters_orf_gene_name(self) -> None:
        """Gene names containing 'orf' are removed."""
        df = self._make_df("1", "protein_coding", "C1orf123")
        result = FetchGencode.filter_protein_coding_genes(df)
        self.assertEqual(len(result), 0)


class TestAddLengthColumn(unittest.TestCase):
    """Test cases for FetchGencode.add_length_column."""

    def test_length_computed(self) -> None:
        """Length column is computed as end - start."""
        df = pd.DataFrame({"start": ["100", "200"], "end": ["300", "500"]})
        result = FetchGencode.add_length_column(df)
        self.assertEqual(result.length.tolist(), [200, 300])

    def test_types_cast_to_int(self) -> None:
        """Start and end columns are cast to int."""
        df = pd.DataFrame({"start": ["100"], "end": ["200"]})
        result = FetchGencode.add_length_column(df)
        self.assertTrue(pd.api.types.is_integer_dtype(result.start))
        self.assertTrue(pd.api.types.is_integer_dtype(result.end))

    def test_does_not_modify_original(self) -> None:
        """Original DataFrame is not mutated."""
        df = pd.DataFrame({"start": ["100"], "end": ["200"]})
        FetchGencode.add_length_column(df)
        self.assertNotIn("length", df.columns)


class TestProcessSingleGene(unittest.TestCase):
    """Test cases for FetchGencode.process_single_gene."""

    def _make_gene_features(self) -> pd.DataFrame:
        """Create a minimal multi-feature gene fixture."""
        gene_id = "ENSG00000001"
        transcript_id = "ENST00000001"
        return pd.DataFrame(
            [
                {
                    "chr": "chr1",
                    "type": "gene",
                    "start": 100,
                    "end": 400,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": "",
                    "transcript_type": "",
                    "length": 300,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
                {
                    "chr": "chr1",
                    "type": "transcript",
                    "start": 100,
                    "end": 400,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": transcript_id,
                    "transcript_type": "protein_coding",
                    "length": 300,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
                {
                    "chr": "chr1",
                    "type": "exon",
                    "start": 100,
                    "end": 200,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": transcript_id,
                    "transcript_type": "protein_coding",
                    "length": 100,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
                {
                    "chr": "chr1",
                    "type": "exon",
                    "start": 300,
                    "end": 400,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": transcript_id,
                    "transcript_type": "protein_coding",
                    "length": 100,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
                {
                    "chr": "chr1",
                    "type": "CDS",
                    "start": 120,
                    "end": 200,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": transcript_id,
                    "transcript_type": "protein_coding",
                    "length": 80,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
                {
                    "chr": "chr1",
                    "type": "CDS",
                    "start": 300,
                    "end": 380,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": transcript_id,
                    "transcript_type": "protein_coding",
                    "length": 80,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
                {
                    "chr": "chr1",
                    "type": "UTR",
                    "start": 100,
                    "end": 120,
                    "strand": "+",
                    "gene_id": gene_id,
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": transcript_id,
                    "transcript_type": "protein_coding",
                    "length": 20,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
            ]
        )

    def test_returns_tuple(self) -> None:
        """Returns a (gene_df, arrow_df) tuple."""
        features = self._make_gene_features()
        result = FetchGencode.process_single_gene("ENSG00000001", features)
        assert result is not None
        gene_df, arrow_df = result
        self.assertIsInstance(gene_df, pd.DataFrame)
        self.assertIsInstance(arrow_df, pd.DataFrame)

    def test_gene_df_has_exon_intron(self) -> None:
        """Gene DataFrame contains both exon and intron features."""
        features = self._make_gene_features()
        result = FetchGencode.process_single_gene("ENSG00000001", features)
        assert result is not None
        gene_df, _ = result
        types = gene_df.type.unique().tolist()
        self.assertIn("exon", types)
        self.assertIn("intron", types)

    def test_arrow_df_has_cds_utr(self) -> None:
        """Arrow DataFrame contains CDS and UTR features."""
        features = self._make_gene_features()
        result = FetchGencode.process_single_gene("ENSG00000001", features)
        assert result is not None
        _, arrow_df = result
        types = arrow_df.type.unique().tolist()
        self.assertIn("CDS", types)
        self.assertIn("UTR", types)

    def test_returns_none_for_no_protein_coding_transcript(self) -> None:
        """Returns None when no protein-coding transcript exists."""
        features = pd.DataFrame(
            [
                {
                    "chr": "chr1",
                    "type": "transcript",
                    "start": 100,
                    "end": 400,
                    "strand": "+",
                    "gene_id": "ENSG00000001",
                    "gene_name": "TESTGENE",
                    "gene_type": "protein_coding",
                    "transcript_id": "ENST00000001",
                    "transcript_type": "lncRNA",
                    "length": 300,
                    "ccdsid": float("nan"),
                    "havana_transcript": float("nan"),
                },
            ]
        )
        result = FetchGencode.process_single_gene("ENSG00000001", features)
        self.assertIsNone(result)

    def test_gene_df_columns(self) -> None:
        """Gene DataFrame contains expected columns."""
        features = self._make_gene_features()
        result = FetchGencode.process_single_gene("ENSG00000001", features)
        assert result is not None
        gene_df, _ = result
        for col in ["chr", "start", "end", "type", "gene_id", "gene_name", "transcript_id"]:
            self.assertIn(col, gene_df.columns)


class TestSampleGtfPipeline(unittest.TestCase):
    """Integration tests running the sample GTF through the full processing pipeline."""

    gencode_raw: pd.DataFrame
    parsed: pd.DataFrame
    stripped: pd.DataFrame
    filtered: pd.DataFrame
    with_length: pd.DataFrame

    @classmethod
    def setUpClass(cls) -> None:
        """Load sample GTF once and run through the full pipeline."""
        # Load raw GTF exactly like retrieve_data does (comment='#' skips header lines):
        raw_tsv = pd.read_csv(
            SAMPLE_GTF, sep="\t", comment="#", header=None, compression="gzip"
        )
        cls.gencode_raw = FetchGencode.parse_raw_gtf(raw_tsv)
        cls.parsed = FetchGencode.parse_gtf_annotations(cls.gencode_raw)
        cls.stripped = FetchGencode.strip_gene_id_version(cls.parsed)
        cls.filtered = FetchGencode.filter_protein_coding_genes(cls.stripped)
        cls.with_length = FetchGencode.add_length_column(cls.filtered)

    # -- parse_raw_gtf ---------------------------------------------------

    def test_raw_shape(self) -> None:
        """Raw GTF has 64 rows and the 6 expected columns."""
        self.assertEqual(len(self.gencode_raw), 64)
        self.assertEqual(
            list(self.gencode_raw.columns),
            ["chr", "type", "start", "end", "strand", "annotation"],
        )

    def test_raw_chromosomes(self) -> None:
        """All rows are on chr8."""
        self.assertEqual(self.gencode_raw.chr.unique().tolist(), ["chr8"])

    def test_raw_feature_types(self) -> None:
        """Raw data contains expected GTF feature types."""
        types = set(self.gencode_raw.type.unique())
        for expected in ("gene", "transcript", "exon", "CDS", "UTR"):
            self.assertIn(expected, types)

    # -- parse_gtf_annotations -------------------------------------------

    def test_parsed_annotation_columns(self) -> None:
        """Parsed DataFrame has gene_id, gene_name, gene_type columns."""
        for col in ("gene_id", "gene_name", "gene_type"):
            self.assertIn(col, self.parsed.columns)

    def test_parsed_annotation_dropped(self) -> None:
        """Raw annotation column is removed after parsing."""
        self.assertNotIn("annotation", self.parsed.columns)

    def test_parsed_gene_name(self) -> None:
        """All rows belong to gene OPRK1."""
        self.assertEqual(self.parsed.gene_name.unique().tolist(), ["OPRK1"])

    def test_parsed_row_count_unchanged(self) -> None:
        """Row count is preserved after annotation parsing."""
        self.assertEqual(len(self.parsed), 64)

    # -- strip_gene_id_version -------------------------------------------

    def test_stripped_gene_id(self) -> None:
        """Version suffix is removed from gene_id."""
        self.assertEqual(
            self.stripped.gene_id.unique().tolist(), ["ENSG00000082556"]
        )

    # -- filter_protein_coding_genes -------------------------------------

    def test_filter_keeps_conventional_chr(self) -> None:
        """Rows on chr8 (conventional chromosome) are kept after filtering."""
        self.assertGreater(len(self.filtered), 0)

    def test_filter_only_protein_coding(self) -> None:
        """All rows after filtering have gene_type protein_coding."""
        self.assertTrue((self.filtered.gene_type == "protein_coding").all())

    def test_filter_single_gene(self) -> None:
        """Only the one gene from the sample survives the filter."""
        gene_count = len(self.filtered.loc[self.filtered.type == "gene"])
        self.assertEqual(gene_count, 1)

    # -- add_length_column -----------------------------------------------

    def test_add_length_column_present(self) -> None:
        """length column is added to the filtered data."""
        self.assertIn("length", self.with_length.columns)

    def test_add_length_values(self) -> None:
        """length equals end - start for every row."""
        self.assertTrue(
            (self.with_length.length == self.with_length.end - self.with_length.start).all()
        )

    # -- process_single_gene ---------------------------------------------

    def test_process_single_gene_from_sample(self) -> None:
        """process_single_gene produces exon/intron and arrow data for OPRK1."""
        groups = list(self.with_length.groupby(["gene_id"]))
        self.assertEqual(len(groups), 1)

        (gene_id,), features = groups[0]
        self.assertEqual(gene_id, "ENSG00000082556")

        result = FetchGencode.process_single_gene(gene_id, features)
        assert result is not None
        gene_df, arrow_df = result

        # Gene DataFrame should have exon + intron features:
        gene_types = set(gene_df.type.unique())
        self.assertIn("exon", gene_types)
        self.assertIn("intron", gene_types)

        # Arrow DataFrame should have CDS + UTR features:
        arrow_types = set(arrow_df.type.unique())
        self.assertIn("CDS", arrow_types)
        self.assertIn("UTR", arrow_types)

        # All output rows should reference OPRK1:
        self.assertTrue((gene_df.gene_name == "OPRK1").all())
        self.assertTrue((arrow_df.gene_name == "OPRK1").all())

        # Gene id should be consistent:
        self.assertTrue((gene_df.gene_id == "ENSG00000082556").all())
        self.assertTrue((arrow_df.gene_id == "ENSG00000082556").all())

    def test_canonical_transcript_selection(self) -> None:
        """Canonical transcript for OPRK1 is selected from CCDS transcripts.

        OPRK1 has 3 protein_coding transcripts with CCDS IDs. The one with
        the longest CDS should be selected as canonical.
        """
        (gene_id,), features = list(self.with_length.groupby(["gene_id"]))[0]
        result = FetchGencode.process_single_gene(gene_id, features)
        assert result is not None
        gene_df, _ = result

        # The canonical transcript should be one of the known OPRK1 transcripts:
        transcript_id = gene_df.transcript_id.unique().tolist()
        self.assertEqual(len(transcript_id), 1)
        self.assertIn(
            transcript_id[0],
            [
                "ENST00000265572.8",
                "ENST00000673285.2",
                "ENST00000524278.5",
                "ENST00000520287.5",
            ],
        )


if __name__ == "__main__":
    unittest.main()
