import os
from typing import List, Optional

import h5py
import numpy as np
from scipy import stats

from artifactsRemoval.Artifact_remove import Artifact_remove


def _save_filtered_hdf5(remover: Artifact_remove, out_path: str, keep_mask: np.ndarray) -> None:
    """Persist the currently kept spots in Visium HDF5 format."""
    matrix = remover.Tissue.tissue_matrix
    if not hasattr(matrix, "tocsc"):
        raise ValueError("tissue_matrix must be a scipy sparse matrix")
    matrix = matrix.tocsc()[:, keep_mask]

    barcodes_series = remover.Tissue.barcode_list.get_column("barcode")
    barcodes = barcodes_series.to_numpy()[keep_mask].astype("S")

    features = remover.Tissue.feature_list
    feature_id = features.get_column("feature_id").to_numpy().astype("S")
    feature_name = features.get_column("feature_name").to_numpy().astype("S")
    if "feature_type" in features.columns:
        feature_type = features.get_column("feature_type").to_numpy().astype("S")
    else:
        feature_type = np.array([b"Gene Expression"] * len(features), dtype="S")
    genome = np.array([b"GRCh38"] * len(features), dtype="S")
    all_tag_keys = np.array([b"feature_type"])

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with h5py.File(out_path, "w") as f:
        grp = f.create_group("matrix")
        grp.create_dataset("barcodes", data=barcodes)
        grp.create_dataset("data", data=matrix.data)
        grp.create_dataset("indices", data=matrix.indices)
        grp.create_dataset("indptr", data=matrix.indptr)
        grp.create_dataset("shape", data=np.array(matrix.shape, dtype=np.int32))

        feat_grp = grp.create_group("features")
        feat_grp.create_dataset("id", data=feature_id)
        feat_grp.create_dataset("name", data=feature_name)
        feat_grp.create_dataset("feature_type", data=feature_type)
        feat_grp.create_dataset("genome", data=genome)
        feat_grp.create_dataset("_all_tag_keys", data=all_tag_keys)


def run_iterative_edge_cleanup(
    data_dir: str,
    output_dir: str,
    p_threshold: float = 0.05,
    initial_edge_distance: int = 1,
    max_iterations: Optional[int] = None,
    remove_border_first: bool = True,
    save_prefix: str = "filtered_feature_bc_matrix",
) -> List[dict]:
    """Iteratively remove border/edge spots until the edge quality matches interior tissue."""
    remover = Artifact_remove(data_dir)
    os.makedirs(output_dir, exist_ok=True)

    barcodes = remover.Tissue.barcode_list.get_column("barcode").to_list()
    depth_df = remover.tissue_depth_df
    depth_map = dict(zip(depth_df["barcode"].to_list(), depth_df["spot_depth"].to_list()))
    depth_values = np.array([depth_map.get(b, 0) for b in barcodes], dtype=int)
    gene_counts = np.asarray(remover.Tissue.tissue_matrix.sum(axis=0)).ravel()

    def inclusion_mask() -> np.ndarray:
        inclusion_map = dict(
            zip(
                remover.spot_inclusion_condition["barcode"].to_list(),
                remover.spot_inclusion_condition["in_tissue"].to_list(),
            )
        )
        return np.array([bool(inclusion_map.get(b, 0)) for b in barcodes], dtype=bool)

    def save_snapshot(label: str) -> None:
        mask = inclusion_mask()
        out_path = os.path.join(output_dir, f"{save_prefix}.{label}.h5")
        _save_filtered_hdf5(remover, out_path, mask)

    if remove_border_first:
        remover.remove_border()
        save_snapshot("border_removed")

    history: List[dict] = []
    max_depth = int(depth_values.max())
    current_distance = max(1, int(initial_edge_distance))
    iteration = 0

    while current_distance <= max_depth:
        mask = inclusion_mask()
        edge_mask = (depth_values == current_distance) & mask
        interior_mask = (depth_values > current_distance) & mask
        if not edge_mask.any() or not interior_mask.any():
            break

        p_val = stats.ttest_ind(
            gene_counts[edge_mask], gene_counts[interior_mask], equal_var=False, nan_policy="omit"
        ).pvalue
        history.append({"edge_distance": int(current_distance), "p_value": float(p_val)})

        if p_val < p_threshold:
            remover.remove_edge(distance=current_distance)
            iteration += 1
            save_snapshot(f"iter{iteration}")
            current_distance += 1
            if max_iterations is not None and iteration >= max_iterations:
                break
        else:
            break

    return history
