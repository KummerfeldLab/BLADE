# BLADE

A lightweight toolkit for detecting and removing tissue artifacts in spatial transcriptomics data.

This repository provides two Python workflows:

- `Artifact_remove` – optimized for standard Visium data (matrix/tsv inputs).
- `Artifact_remove_HD` – optimized for Visium HD data (HDF5 + Parquet/CSV inputs).

Both workflows share a similar API: border removal, edge removal, malfunction detection, and utilities for reviewing or saving masks.

## Installation

Clone the repository and install it locally (editable mode recommended during development):

```
git clone https://github.com/KummerfeldLab/BLADE.git
cd BLADE
pip install -e .
```

Dependencies are listed in `setup.py` (NumPy, pandas, SciPy, matplotlib, h5py, pyarrow for HD Parquet support).

## Quick start

```python
from artifactsRemoval import Artifact_remove, Artifact_remove_HD

# Standard Visium data
visium = Artifact_remove(dir="/path/to/visium/sample")
visium.remove_border()
visium.remove_edge(distance=2)
visium.remove_malfunction(use_clusters=False)
visium.save("/path/to/output")

# Visium HD data
hd = Artifact_remove_HD(
    data_dir="/path/to/Visium_HD/sample",
    neighbor_mode="k",
    k_neighbors=8,
)
hd.remove_border(distance=float(hd._grid_radius) if hd._grid_radius else 0.0)
hd.remove_edge(distance=2)
hd.remove_malfunction()  # clustered removal by default
hd.save("/path/to/output")
```

### Choosing the right class

| Data type | Class to use | Notes |
| --- | --- | --- |
| Visium (matrix/tsv) | `Artifact_remove` | Uses lattice-based neighbor lookups. |
| Visium HD (HDF5) | `Artifact_remove_HD` | Uses KD-tree spatial neighbors and supports Parquet/CSV positions. Requires `pyarrow` for Parquet. |

## Tutorials

- **HD workflow notebook** – `notebooks/Visium_HD_Artifact_Cleanup.ipynb` demonstrates the HD pipeline end-to-end.
- Legacy tutorials for standard Visium remain under `artifactsRemoval/tutorial_package.py` and the `old_ver/` directory.

## Packaging notes

- `setup.py` exposes the package as `artifactsRemoval`. Install with `pip install -e .` during development or `pip install git+https://github.com/KummerfeldLab/BLADE.git` for users.
- `artifactsRemoval/__init__.py` exports both `Artifact_remove` and `Artifact_remove_HD`.

## Development workflow

1. Create a feature branch from `lw_tissue_composition` (as done for this integration).
2. Implement changes (e.g., new HD functionality).
3. Run the notebooks/tests in `notebooks/` or `artifactsRemoval/test_*` to verify behavior.
4. Commit and push your branch, then open a pull request.

## Support and limitations

- Artifact definitions assume Visium-style grids; adjust neighbor radius and cluster thresholds as needed for other technologies.
- Cluster-aware malfunction removal (`remove_malfunction`) removes only contiguous outlier patches with size ≥ 15 by default; set `use_clusters=False` to revert to the original per-spot removal.
- Ensure `pyarrow` is installed when working with Parquet position files for Visium HD.
