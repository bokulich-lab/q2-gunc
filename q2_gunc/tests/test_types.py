# ----------------------------------------------------------------------------
# Copyright (c) 2025, Bokulich Lab.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os

import pandas as pd
import pytest
from q2_types.reference_db import ReferenceDB
from qiime2.core.exceptions import ValidationError
from qiime2.plugin.testing import TestPluginBase
import qiime2

from q2_gunc.types import (
    GUNCResultsFormat,
    GUNCGeneCountsFormat,
    GUNCHTMLPlotFormat,
    GUNCResultsDirectoryFormat,
    GUNCDatabaseDirFmt,
    GUNCDB,
    GUNCResults,
)


class TestTypes(TestPluginBase):
    package = "q2_gunc.tests"

    def test_guncresultsformat_valid(self):
        f = self.get_data_path("valid_gunc.tsv")
        fmt = GUNCResultsFormat(f, mode="r")
        fmt.validate()

    def test_guncresultsformat_invalid(self):
        f = self.get_data_path("invalid_gunc.tsv")
        fmt = GUNCResultsFormat(f, mode="r")
        with pytest.raises(ValidationError):
            fmt.validate()

    def test_guncgenecountsformat_valid(self):
        f = self.get_data_path("valid_gene_counts.json")
        fmt = GUNCGeneCountsFormat(f, mode="r")
        fmt.validate()

    def test_guncgenecountsformat_invalid(self):
        f = self.get_data_path("invalid_gene_counts.json")
        fmt = GUNCGeneCountsFormat(f, mode="r")
        with pytest.raises(ValidationError):
            fmt.validate()

    def test_gunchtmlplotformat_valid(self):
        f = self.get_data_path("valid_plot.html")
        fmt = GUNCHTMLPlotFormat(f, mode="r")
        fmt.validate()

    def test_gunchtmlplotformat_invalid(self):
        f = self.get_data_path("invalid_plot.html")
        fmt = GUNCHTMLPlotFormat(f, mode="r")
        # BeautifulSoup is lenient, so this may not error, but we check anyway
        try:
            fmt.validate()
        except ValidationError:
            pass

    def test_guncresultsdirectoryformat_no_samples(self):
        f = self.get_data_path("results")
        fmt = GUNCResultsDirectoryFormat(f, mode="r")
        fmt.validate()

        obs_file_dict = fmt.file_dict()
        exp_file_dict = {"": str(f)}
        self.assertDictEqual(obs_file_dict, exp_file_dict)

    def test_guncresultsdirectoryformat_with_samples(self):
        f = self.get_data_path("results-per-sample")
        fmt = GUNCResultsDirectoryFormat(f, mode="r")
        fmt.validate()

        obs_file_dict = fmt.file_dict()
        exp_file_dict = {
            "SRR9640343": os.path.join(f, "SRR9640343"),
            "SRR9640344": os.path.join(f, "SRR9640344"),
        }
        self.assertDictEqual(obs_file_dict, exp_file_dict)

    def test_guncresultsdirectoryformat_no_plots(self):
        f = self.get_data_path("results-no-plots")
        fmt = GUNCResultsDirectoryFormat(f, mode="r")
        fmt.validate()

    def test_guncdatabasedirfmt(self):
        f = self.get_data_path("db")
        fmt = GUNCDatabaseDirFmt(f, mode="r")
        fmt.validate()

    def test_db_semantic_type_registration(self):
        self.assertRegisteredSemanticType(GUNCDB)

    def test_results_semantic_type_registration(self):
        self.assertRegisteredSemanticType(GUNCResults)

    def test_db_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            ReferenceDB[GUNCDB], GUNCDatabaseDirFmt
        )

    def test_results_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            GUNCResults, GUNCResultsDirectoryFormat
        )

    def test_transformer_registration_dataframe(self):
        """Test that DataFrame transformer is registered."""
        import pandas as pd
        self.assertTransformerRegistered(GUNCResultsDirectoryFormat, pd.DataFrame)

    def test_transformer_registration_metadata(self):
        """Test that Metadata transformer is registered."""
        from qiime2 import Metadata
        self.assertTransformerRegistered(GUNCResultsDirectoryFormat, Metadata)

    def test_transformer_no_samples(self):
        """Test transformer with non-partitioned data (FeatureData[MAG])."""
        f = self.get_data_path("results")
        fmt = GUNCResultsDirectoryFormat(f, mode="r")
        
        # Test DataFrame transformer
        df = fmt.view(view_type=pd.DataFrame)
        
        # Should have composite indices with name "id"
        self.assertEqual(df.index.name, "id")
        
        # Should contain composite indices like "genome_id_taxonomic_level_index"
        # Each MAG should have 7 rows (one for each taxonomic level)
        expected_mags = ["0c20367d-4775-43f1-90c6-1a36afc5e4da", "1da59757-769b-4713-923d-e3d2e60690c9"]
        genome_ids_in_index = set([idx.split("_")[0] for idx in df.index])
        self.assertEqual(genome_ids_in_index, set(expected_mags))
        
        # Should have 14 rows total (2 MAGs × 7 taxonomic levels)
        self.assertEqual(len(df), 14)
        
        # Should not have a sample_id column for non-partitioned data
        self.assertNotIn("sample_id", df.columns)
        
        # Should have the GUNC result columns
        expected_columns = [
            "genome", "n_genes_called", "n_genes_mapped", "n_contigs", "taxonomic_level",
            "proportion_genes_retained_in_major_clades", "genes_retained_index",
            "clade_separation_score", "contamination_portion", "n_effective_surplus_clades",
            "mean_hit_identity", "reference_representation_score", "pass.GUNC"
        ]
        for col in expected_columns:
            self.assertIn(col, df.columns)

        # Test Metadata transformer
        metadata = fmt.view(view_type=qiime2.Metadata)
        self.assertIsInstance(metadata, qiime2.Metadata)
        
        # Should have same structure as DataFrame
        metadata_df = metadata.to_dataframe()
        self.assertEqual(metadata_df.index.name, "id")
        self.assertEqual(len(metadata_df), 14)

    def test_transformer_with_samples(self):
        """Test transformer with partitioned data (SampleData[MAGs])."""
        f = self.get_data_path("results-per-sample")
        fmt = GUNCResultsDirectoryFormat(f, mode="r")
        
        # Test DataFrame transformer
        df = fmt.view(view_type=pd.DataFrame)
        
        # Should have composite indices with name "id"
        self.assertEqual(df.index.name, "id")
        
        # Should have a sample_id column for partitioned data
        self.assertIn("sample_id", df.columns)
        
        # Should contain MAGs from both samples
        expected_samples = ["SRR9640343", "SRR9640344"]
        observed_samples = df["sample_id"].unique().tolist()
        self.assertEqual(sorted(observed_samples), sorted(expected_samples))
        
        # Should have 4 MAGs × 7 taxonomic levels = 28 rows
        self.assertEqual(len(df), 28)
        
        # Should have the GUNC result columns plus sample_id
        expected_columns = [
            "genome", "n_genes_called", "n_genes_mapped", "n_contigs", "taxonomic_level",
            "proportion_genes_retained_in_major_clades", "genes_retained_index",
            "clade_separation_score", "contamination_portion", "n_effective_surplus_clades",
            "mean_hit_identity", "reference_representation_score", "pass.GUNC", "sample_id"
        ]
        for col in expected_columns:
            self.assertIn(col, df.columns)

        # Test Metadata transformer
        metadata = fmt.view(view_type=qiime2.Metadata)
        self.assertIsInstance(metadata, qiime2.Metadata)
        
        # Should have same structure as DataFrame
        metadata_df = metadata.to_dataframe()
        self.assertEqual(metadata_df.index.name, "id")
        self.assertEqual(len(metadata_df), 28)
        self.assertIn("sample_id", metadata_df.columns)
