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
A lightweight toolkit for detecting and removing tissue artifacts in spatial transcriptomics data.

This repository provides two Python workflows:

- `Artifact_remove` – optimized for standard Visium data (matrix/tsv inputs).
- `Artifact_remove_HD` – optimized for Visium HD data (HDF5 + Parquet/CSV inputs).

Both workflows share a similar API: border removal, edge removal, malfunction detection, and utilities for reviewing or saving masks.

## Installation

Install directly from GitHub:

```
pip install git+https://github.com/KummerfeldLab/BLADE.git
```

For editable installs (development):

```
git clone https://github.com/KummerfeldLab/BLADE.git
cd BLADE
pip install -e .
```

Dependencies include NumPy, pandas, SciPy, matplotlib, h5py, and pyarrow (for Parquet support in HD workflows).

## Usage examples

### Standard Visium (matrix/tsv inputs)

```python
from artifactsRemoval import Artifact_remove

dir_path = "/path/to/visium/sample"
visium = Artifact_remove(dir=dir_path)

visium.remove_border()
visium.remove_edge(distance=2)
visium.remove_malfunction(use_clusters=False)

visium.review_removing()
visium.save("/path/to/output")
```

### Visium HD (HDF5 + Parquet/CSV)

```python
from artifactsRemoval import Artifact_remove_HD

hd = Artifact_remove_HD(
    data_dir="/path/to/Visium_HD/sample",
    positions_path="/path/to/Visium_HD/sample/spatial/tissue_positions.parquet",
    h5_path="/path/to/Visium_HD/sample/filtered_feature_cell_matrix.h5",
    neighbor_mode="k",
    k_neighbors=8,
)

border_width = float(hd._grid_radius) if hd._grid_radius else 0.0
hd.remove_border(distance=border_width)
hd.remove_edge(distance=2)
hd.remove_malfunction()  # cluster-aware by default (min cluster size 15)

hd.save("/path/to/output")
hd.save_hdf5("/path/to/output/filtered_feature_cell_matrix.filtered.h5")
```



## One-button edge cleanup

To automatically iterate border/edge removal until edge quality matches the interior, call `run_iterative_edge_cleanup`. The helper loads the sample, optionally removes the border, runs two-sided t-tests comparing the current edge to interior counts, removes one edge layer when `p < threshold`, and saves an HDF5 snapshot after each removal.

```python
from artifactsRemoval import run_iterative_edge_cleanup

history = run_iterative_edge_cleanup(
    data_dir="/path/to/visium/sample",
    output_dir="/tmp/edge_runs",
    p_threshold=0.05,
    initial_edge_distance=1,
    max_iterations=5,
    remove_border_first=True,
)
print(history)
# Saves files: filtered_feature_bc_matrix.border_removed.h5, filtered_feature_bc_matrix.iter1.h5, ...
```

### Usage example (iterative deletion loop)

```
from artifactsRemoval import run_iterative_edge_cleanup

# Automatically peel edge layers until the edge no longer differs from the interior
history = run_iterative_edge_cleanup(
    data_dir="/data/spaceranger/sample1",
    output_dir="/data/results/sample1_cleanup",
    p_threshold=0.05,          # stop when edge vs interior p-value is >= 0.05
    initial_edge_distance=1,   # start with distance-1 edge layer
    max_iterations=None,       # remove until convergence or tissue is exhausted
    remove_border_first=True,  # drop annotated border before testing edges
)

for step in history:
    print(f"Edge distance {step['edge_distance']} p={step['p_value']:.3g}")
# HDF5 snapshots land in output_dir after each removal (iterN.h5) plus the optional border_removed.h5
```

## Tutorials

- `notebooks/Visium_HD_Artifact_Cleanup.ipynb` – interactive walkthrough for the Visium HD workflow.
- `artifactsRemoval/tutorial_package.py` – script illustrating standard Visium usage with `Artifact_remove`.

## Packaging notes

- The package installs as `artifactsRemoval` and exports both `Artifact_remove` and `Artifact_remove_HD`.
- Ensure `pyarrow` is installed when working with Parquet position files in Visium HD data.
