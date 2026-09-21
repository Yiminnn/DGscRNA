"""Regression tests for release review findings; use SLURM at the project site."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from dgscrna.reference.io import marker_coverage
from dgscrna.reference.runner import run_reference


class ReferenceContractTests(unittest.TestCase):
    def test_zero_overlap_cannot_create_single_panel_seed_calls(self):
        with self.assertRaisesRegex(ValueError, "zero overlap"):
            marker_coverage({"brain": {"Astrocyte": ["GFAP"]}}, ["ENSG00001"])

    def test_unselected_library_is_reported_without_blocking_requested_library(self):
        result = marker_coverage({"brain": {"Astrocyte": ["GFAP", "GFAP", "AQP4"]},
                                  "other": {"Absent": ["ZZZ"]}}, ["GFAP"], "brain")
        self.assertEqual(result["brain"]["panels"]["Astrocyte"], {"denominator": 3, "retained_unique": 1})
        self.assertEqual(result["other"]["retained_marker_genes"], 0)

    def test_api_seed_cannot_disagree_with_r_path_conversion(self):
        for seed in ["42", 42.5, True, -1, 2**31]:
            with self.subTest(seed=seed), self.assertRaisesRegex(ValueError, "seed"):
                run_reference(counts="unused", markers="unused", out="unused", seed=seed)

    def test_export_preserves_numeric_and_na_cell_identifiers(self):
        from dgscrna.reference.runner import _export
        import pandas as pd
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            route = root / "prep/route"
            route.mkdir(parents=True)
            cells = ["001", "002", "NA"]
            pd.DataFrame({"cell_id": cells, "final090": ["A", "A", "B"]}).to_csv(root / "predictions.csv.gz", index=False)
            pd.DataFrame({"cell_id": cells, "cluster": ["0", "0", "1"]}).to_csv(route / "clusters.csv", index=False)
            pd.DataFrame({"UMAP_1": [1., 2., 3.], "UMAP_2": [4., 5., 6.]}, index=cells).to_csv(root / "prep/UMAP2.csv")
            record = dict(predictions="predictions.csv.gz", route="route", arm="L00_mean", library="brain",
                          cutoff="mean", dl_status="no_op_all_initially_known", training_executed=False)
            _export(root, root / "prep", [record])
            for name in ["annotations.csv.gz", "embedding.csv"]:
                frame = pd.read_csv(root / name, dtype=str, keep_default_na=False)
                self.assertEqual(frame.cell_id.tolist(), cells)
