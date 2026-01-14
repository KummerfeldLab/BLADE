# BLADE (R)

This is the R implementation of BLADE, built on top of the **polars** R bindings.
It mirrors the Python workflow: identify border/edge spots using neighbor graphs,
remove artifacts, and iterate edge cleanup until edge quality matches the interior.

## Structure

- `R/` contains package source files.
- `man/` contains documentation stubs for CRAN.

## Quick start

```r
library(BLADE)

tissue <- blade_read_visium(
  "/path/to/visium/sample",
  h5_path = "/path/to/visium/sample/filtered_feature_bc_matrix/filtered_feature_bc_matrix.h5"
)

tissue <- blade_remove_border(tissue)

tissue <- blade_remove_edge(tissue, distance = 2)

tissue <- blade_remove_malfunction(tissue)

blade_save_positions(tissue, "/tmp/tissue_positions.cleaned.csv")
```

## Visium HD iterative cleanup

```r
history <- blade_run_iterative_edge_cleanup_hd(
  data_dir = "/path/to/Visium_HD/sample",
  output_dir = "/tmp/hd_edge_runs",
  p_threshold = 0.05,
  initial_edge_distance = 1,
  max_iterations = 5,
  remove_border_first = TRUE
)
print(history)
```

## HDF5 input notes

- `blade_read_visium()` prefers `filtered_feature_bc_matrix.h5` when available.
- If the HDF5 file is missing, it falls back to the MTX/TSV inputs.

## Performance notes (runtime estimates)

The neighbor graph is built with a grid-hash approach, so runtime is dominated by
local candidate checks rather than an all-pairs scan.

- Neighbor indexing: roughly **O(n * m)**, where *m* is the average number of
  candidate points in nearby grid cells. For typical Visium layouts, *m* is small
  (tens of points), giving near-linear scaling.
- Edge/border depth propagation: **O(n + e)**, where *e* is the number of neighbor
  links.
- Iterative edge cleanup: **O(k * (n + e))**, where *k* is the number of edge layers
  removed.

If runtime becomes a bottleneck, the fastest improvements are:

1. Precompute and cache neighbor maps at the target radius.
2. Use a smaller neighbor radius when appropriate (reduces candidate checks).
3. Run on a subset of barcodes if only coarse edge detection is needed.
