"""Boundary tests for the public reference adapters; execute inside SLURM.

These tests reconstruct the exported R sparse input and distinguish raw counts
from misleading alternate inputs. They do not fit a scientific model.
"""
from pathlib import Path
import gzip
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Reference numerical adapter tests must run inside SLURM")

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from dgscrna.reference.io import counts_files, identifier, load_markers, prepare_counts, sha
from dgscrna.reference.runner import _checked_receipt


class ReferenceInputTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        # Distinct cell/gene values expose a transpose, reorder or wrong ID column.
        self.raw = (1 + np.arange(40 * 36).reshape(40, 36) % 17).astype(float)
        self.raw[2:, 0] = 0  # This gene has only two detected cells and must disappear.
        self.raw[3:, 1] = 0  # Exactly three detected cells must retain this gene.
        self.cells = [f"cell-{i:02d}" for i in range(40)]
        self.genes = [f"SYMBOL{i:02d}" for i in range(36)]

    def tenx(self, name="counts", *, raw=None, compressed=False, three_columns=False,
             genes=None, cells=None, assay="Gene Expression"):
        directory = self.root / name
        directory.mkdir()
        x = sparse.coo_matrix(self.raw if raw is None else raw)
        matrix = directory / "matrix.mtx"
        mmwrite(matrix, x.T, field="real")
        features = [f"ENSG{i:06d}\t{gene}" + (f"\t{assay}" if three_columns else "")
                    for i, gene in enumerate(self.genes if genes is None else genes)]
        gene_file = directory / ("features.tsv" if three_columns else "genes.tsv")
        gene_file.write_text("\n".join(features) + "\n")
        barcodes = directory / "barcodes.tsv"
        barcodes.write_text("\n".join(self.cells if cells is None else cells) + "\n")
        if compressed:
            for path in [matrix, gene_file, barcodes]:
                with gzip.open(str(path) + ".gz", "wb") as stream:
                    stream.write(path.read_bytes())
                path.unlink()
        return directory

    def exported(self, directory, **kwargs):
        output = self.root / (directory.stem + "_export")
        manifest = prepare_counts(directory, output, "example", **kwargs)
        matrix = sparse.csr_matrix(
            (np.fromfile(output / "x.bin", dtype="<f8"),
             np.fromfile(output / "i.bin", dtype="<i4"),
             np.fromfile(output / "p.bin", dtype="<i4")),
            shape=(manifest["n_cells"], manifest["n_genes"]),
        )
        return output, manifest, matrix

    def test_two_and_three_column_gzip_inputs_export_identical_r_counts(self):
        first, m1, x1 = self.exported(self.tenx("two"))
        second, m2, x2 = self.exported(self.tenx("three", compressed=True, three_columns=True))
        np.testing.assert_array_equal(x1.toarray(), self.raw[:, 1:])
        np.testing.assert_array_equal(x2.toarray(), self.raw[:, 1:])
        self.assertEqual(m1["n_genes_before_filter"], 36)
        self.assertEqual(m1["n_genes"], 35)
        self.assertEqual(m1["fitting_files"], m2["fitting_files"])
        self.assertEqual(pd.read_csv(first / "genes.csv").gene.tolist(), self.genes[1:])
        self.assertEqual(pd.read_csv(second / "cells_fit.csv").cell_id.tolist(), self.cells)
        self.assertEqual((first / "INPUT_COMPLETE").read_text().strip(), sha(first / "input_manifest.json"))

    def test_multimodal_feature_rows_are_not_silently_mixed_with_rna(self):
        counts = self.tenx(three_columns=True, assay="Antibody Capture")
        with self.assertRaisesRegex(ValueError, "Gene Expression"):
            self.exported(counts)

    def test_ambiguous_alternate_feature_files_are_rejected(self):
        counts = self.tenx()
        (counts / "features.tsv").write_bytes((counts / "genes.tsv").read_bytes())
        with self.assertRaisesRegex(ValueError, "exactly one genes"):
            counts_files(counts)

    def test_counts_layer_cannot_be_ignored_for_matrixmarket(self):
        with self.assertRaisesRegex(ValueError, "only to h5ad"):
            self.exported(self.tenx(), counts_layer="counts")

    def test_negative_fractional_and_nonfinite_counts_are_rejected(self):
        for index, invalid in enumerate([-1, 0.25, np.nan, np.inf]):
            with self.subTest(value=invalid):
                x = self.raw.copy()
                x[5, 6] = invalid
                with self.assertRaisesRegex(ValueError, "finite nonnegative integer counts"):
                    self.exported(self.tenx(f"invalid{index}", raw=x))

    def test_duplicate_barcodes_and_duplicate_gene_symbols_are_rejected(self):
        cells = self.cells.copy()
        cells[-1] = cells[0]
        with self.assertRaisesRegex(ValueError, "cell identifiers"):
            self.exported(self.tenx("dupe_cells", cells=cells))
        genes = self.genes.copy()
        genes[-1] = genes[0]
        with self.assertRaisesRegex(ValueError, "gene identifiers"):
            self.exported(self.tenx("dupe_genes", genes=genes))

    def test_empty_cell_after_detection_filter_is_rejected(self):
        x = self.raw.copy()
        x[0] = 0
        x[0, 0] = 3
        with self.assertRaisesRegex(ValueError, "positive counts after"):
            self.exported(self.tenx(raw=x))

    def test_dimensions_do_not_silently_truncate_identifiers(self):
        with self.assertRaisesRegex(ValueError, "dimensions"):
            self.exported(self.tenx(cells=self.cells[:-1]))

    def test_small_input_fails_before_r_pca(self):
        with self.assertRaisesRegex(ValueError, "PCA30"):
            self.exported(self.tenx(raw=self.raw[:31], cells=self.cells[:31]))

    def test_duplicate_sparse_entries_do_not_hide_invalid_counts(self):
        base = sparse.coo_matrix(self.raw)
        # -0.25 + 0.25 cancels during COO -> CSR conversion. The raw input is
        # nevertheless fractional and must not be admitted as integer counts.
        x = sparse.coo_matrix(
            (np.r_[base.data, -0.25, 0.25],
             (np.r_[base.row, 5, 5], np.r_[base.col, 6, 6])),
            shape=base.shape,
        )
        counts = self.tenx(raw=x)
        with self.assertRaisesRegex(ValueError, "finite nonnegative integer counts"):
            self.exported(counts)

    def test_large_integer_counts_are_not_silently_rounded(self):
        x = self.raw.copy()
        x[5, 6] = 2**24 + 1  # First positive integer not exactly representable as float32.
        counts = self.tenx(raw=x)
        try:
            _, _, exported = self.exported(counts)
        except ValueError as error:
            self.assertRegex(str(error), "exact|precision|represent|float32")
        else:
            # The public integer-count contract permits an explicit rejection,
            # but not an unnoticed change of count values during conversion.
            self.assertEqual(exported[5, 5], x[5, 6])

    @unittest.skipUnless(importlib.util.find_spec("anndata"), "Optional anndata dependency is unavailable")
    def test_h5ad_requires_explicit_raw_layer_and_preserves_it(self):
        import anndata

        adata = anndata.AnnData(
            X=sparse.csr_matrix(self.raw / 2),
            obs=pd.DataFrame(index=self.cells), var=pd.DataFrame(index=self.genes),
        )
        adata.layers["counts"] = sparse.csr_matrix(self.raw)
        # Truth metadata and normalized .X must not replace the chosen count layer.
        adata.obs["author_label"] = "do_not_use"
        path = self.root / "input.h5ad"
        adata.write_h5ad(path)
        with self.assertRaisesRegex(ValueError, "explicitly set"):
            prepare_counts(path, self.root / "missing_choice", "example")
        with self.assertRaisesRegex(ValueError, "Missing count layer"):
            prepare_counts(path, self.root / "missing_layer", "example", "missing")
        with self.assertRaisesRegex(ValueError, "integer counts"):
            prepare_counts(path, self.root / "normalized_X", "example", "X")
        output, manifest, result = self.exported(path, counts_layer="counts")
        np.testing.assert_array_equal(result.toarray(), self.raw[:, 1:])
        self.assertEqual(pd.read_csv(output / "cells_fit.csv").columns.tolist(), ["cell_id", "batch"])
        self.assertTrue(manifest["truth_labels_excluded_from_fit"])


class MarkerAndIdentifierTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def test_single_library_and_duplicate_genes_preserve_panel_denominator(self):
        expected = {"brain": {"Astrocyte": ["GFAP", "GFAP", "AQP4"]}}
        json_path = self.root / "markers.json"
        json_path.write_text(json.dumps(expected))
        tsv_path = self.root / "markers.tsv"
        tsv_path.write_text("library\tcell_type\tgene\nbrain\tAstrocyte\tGFAP\n"
                            "brain\tAstrocyte\tGFAP\nbrain\tAstrocyte\tAQP4\n")
        self.assertEqual(load_markers(json_path), expected)
        self.assertEqual(load_markers(tsv_path), expected)

    def test_reserved_terminal_states_cannot_be_marker_cell_types(self):
        for panel in ["Unknown", "Undecided"]:
            with self.subTest(panel=panel):
                path = self.root / "markers.json"
                path.write_text(json.dumps({"brain": {panel: ["GFAP"]}}))
                with self.assertRaisesRegex(ValueError, "reserved"):
                    load_markers(path)

    def test_duplicate_json_keys_are_rejected_at_both_levels(self):
        for raw in ['{"brain":{"A":["X"]},"brain":{"B":["Y"]}}',
                    '{"brain":{"A":["X"],"A":["Y"]}}']:
            with self.subTest(raw=raw):
                path = self.root / "markers.json"
                path.write_text(raw)
                with self.assertRaisesRegex(ValueError, "Duplicate marker JSON key"):
                    load_markers(path)

    def test_invalid_sample_identifiers_cannot_escape_output_directory(self):
        for name in ["../sample", "/tmp/sample", "a/b", ".", "..", "has space", "", "x\n"]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                identifier(name)
        self.assertEqual(identifier("TKU-3186.1_A"), "TKU-3186.1_A")


class ReferenceReceiptTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.artifact = self.root / "stage" / "features.txt"
        self.artifact.parent.mkdir()
        self.artifact.write_text("AQP4\nGFAP\n")
        self.receipt = self.root / "receipt.json"
        self.receipt.write_text(json.dumps({
            "plan_sha256": "original-plan", "artifacts": {"stage/features.txt": sha(self.artifact)},
        }))

    def test_checkpoint_verifies_content_not_only_file_existence(self):
        self.assertTrue(_checked_receipt(self.receipt, self.root, "original-plan"))
        self.artifact.write_text("AQP4\nMBP\n")
        with self.assertRaisesRegex(ValueError, "missing or modified"):
            _checked_receipt(self.receipt, self.root, "original-plan")

    def test_deleted_artifact_and_changed_plan_cannot_resume(self):
        with self.assertRaisesRegex(ValueError, "configuration changed"):
            _checked_receipt(self.receipt, self.root, "other-plan")
        self.artifact.unlink()
        with self.assertRaisesRegex(ValueError, "missing or modified"):
            _checked_receipt(self.receipt, self.root, "original-plan")
        self.assertFalse(_checked_receipt(self.root / "absent.json", self.root, "original-plan"))


class LightweightEntryPointTests(unittest.TestCase):
    def test_help_and_public_reference_import_do_not_import_scanpy_or_torch(self):
        program = """
import importlib.abc
import sys
class RejectLegacyImports(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'scanpy', 'torch', 'numpy', 'anndata'}:
            raise RuntimeError('Heavy dependency imported during --help: ' + fullname)
sys.meta_path.insert(0, RejectLegacyImports())
from dgscrna import run_reference
from dgscrna.reference.cli import main
try:
    main(['run', '--help'])
except SystemExit as error:
    assert error.code == 0
else:
    raise AssertionError('argparse help did not exit')
"""
        completed = subprocess.run([sys.executable, "-s", "-c", program], text=True, capture_output=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("--counts-layer", completed.stdout)
        self.assertIn("--stop-after", completed.stdout)


if __name__ == "__main__":
    unittest.main()
