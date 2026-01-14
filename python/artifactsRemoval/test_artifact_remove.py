import os
import unittest
import numpy as np
import polars as pl
import scipy.sparse as sp
import h5py
from artifactsRemoval.Artifact_remove import Artifact_remove


class TestArtifactRemove(unittest.TestCase):
    def setUp(self):
        n_genes, n_barcodes = 5, 3
        data = np.array([1, 2, 3, 4, 5, 6])
        indices = np.array([0, 2, 2, 0, 1, 2])
        indptr = np.array([0, 2, 3, 6])
        shape = (n_genes, n_barcodes)
        self.matrix = sp.csc_matrix((data, indices, indptr), shape=shape)
        self.barcodes = pl.DataFrame({"barcode": [f"BC{i}" for i in range(n_barcodes)]})
        self.features = pl.DataFrame(
            {
                "feature_id": [f"GENE{i}" for i in range(n_genes)],
                "feature_name": [f"GeneName{i}" for i in range(n_genes)],
                "feature_type": ["Gene Expression"] * n_genes,
            }
        )

        class DummyArtifactRemove(Artifact_remove):
            def __init__(self, matrix, barcodes, features):
                self.Tissue = type("Dummy", (), {})()
                self.Tissue.tissue_matrix = matrix
                self.Tissue.barcode_list = barcodes
                self.Tissue.feature_list = features

        self.dummy = DummyArtifactRemove(self.matrix, self.barcodes, self.features)

    def test_save_hdf5(self):
        import tempfile

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test_filtered_feature_bc_matrix.h5")
            self.dummy.save_hdf5(out_path)
            with h5py.File(out_path, "r") as f:
                self.assertIn("matrix", f)
                grp = f["matrix"]
                np.testing.assert_array_equal(grp["barcodes"][:], self.barcodes["barcode"].to_numpy().astype("S"))
                np.testing.assert_array_equal(grp["data"][:], self.matrix.data)
                np.testing.assert_array_equal(grp["indices"][:], self.matrix.indices)
                np.testing.assert_array_equal(grp["indptr"][:], self.matrix.indptr)
                np.testing.assert_array_equal(grp["shape"][:], np.array(self.matrix.shape, dtype=np.int32))
                feat_grp = grp["features"]
                np.testing.assert_array_equal(
                    feat_grp["id"][:], self.features["feature_id"].to_numpy().astype("S")
                )
                np.testing.assert_array_equal(
                    feat_grp["name"][:], self.features["feature_name"].to_numpy().astype("S")
                )
                np.testing.assert_array_equal(
                    feat_grp["feature_type"][:], self.features["feature_type"].to_numpy().astype("S")
                )
                np.testing.assert_array_equal(
                    feat_grp["genome"][:], np.array([b"GRCh38"] * self.matrix.shape[0])
                )
                np.testing.assert_array_equal(feat_grp["_all_tag_keys"][:], np.array([b"feature_type"]))


if __name__ == "__main__":
    unittest.main()
