\name{blade_run_iterative_edge_cleanup_hd}
\alias{blade_run_iterative_edge_cleanup_hd}
\title{Iterative edge cleanup for Visium HD}
\usage{
blade_run_iterative_edge_cleanup_hd(data_dir = NULL, h5_path = NULL,
  positions_path = NULL, output_dir, p_threshold = 0.05,
  initial_edge_distance = 1, max_iterations = NULL,
  remove_border_first = TRUE)
}
\arguments{
  \item{data_dir}{Optional data directory for auto-discovery.}
  \item{h5_path}{Path to a 10x HDF5 matrix file.}
  \item{positions_path}{Path to tissue_positions.parquet or CSV.}
  \item{output_dir}{Directory to save tissue_positions snapshots.}
  \item{p_threshold}{Stop when edge vs interior p-value >= threshold.}
  \item{initial_edge_distance}{Starting edge distance.}
  \item{max_iterations}{Optional maximum iterations.}
  \item{remove_border_first}{Whether to remove border before edge testing.}
}
\value{A data.frame with edge distances and p-values.}
\description{Iteratively removes edge layers for Visium HD until edge matches interior tissue.}
