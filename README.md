## **BLADE**: a package for the purpose of tissue artifacts detection (Edge, Border, abnormal spots). 
It has both Python and R implementations. For detailed usage tutorial refer to usage tutorial section. 

### Main use:  Detecting edge artifacts, border artifacts, and location malfunction artifacts

### Context of use: Spatial transcriptomics data, especially Visium data

### Pitfalls: 
- Use of this method does not guarantee that there are no other types of artifacts.

- For different technologies or tissue samples, some parameters may need to be changed by hand to ensure optimal performance.

- Artifacts may not be detected if they: (1) occur only as gradual changes over larger pieces of tissue, (2) are low in strength, (3) occur in only small portions of the border or tissue edge.

- BLADE has not been thoroughly evaluated for data from technologies other than Visium.

### Principle of operation 

- Edge artifacts and border artifacts are detected by first identifying edge spots and border spots based on the image, and then t-tests are used to compare those spots to interior spots.

- Location malfunctions are identified at the batch level by identifying outlier spots in each image in the batch, identifying spot locations that are outliers in most or all of the images in the batch, and then testing whether there is a connected group of locations that is larger than should occur by chance.


### Theoretical properties and empirical evidence
- Edge and border artifacts detection inherits the theoretical properties from the t-test used.

- Location artifact detection inherits the theoretical properties of the outlier spot detection method. The test for detecting if the largest outlier cluster is implausibly large relies on simulating random spatial distributions tailored to the shape of the tissue samples in the images, and is asymptotically correct as the number of simulated tissue samples increases.

### Best practices

- Images tested with BLADE should be investigated in other ways (e.g. heatmaps, plotting sequencing depth by edge distance, etc.) to assess if any of the pitfall conditions might apply.

- BLADE should be used on Visium and other spatial transcriptomics data prior to other forms of analysis, to prevent bias
----------------------

## Intro
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
