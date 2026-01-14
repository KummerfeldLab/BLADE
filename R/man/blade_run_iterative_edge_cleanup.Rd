\name{blade_run_iterative_edge_cleanup}
\alias{blade_run_iterative_edge_cleanup}
\title{Iterative edge cleanup}
\usage{
blade_run_iterative_edge_cleanup(data_dir, output_dir, p_threshold = 0.05,
  initial_edge_distance = 1, max_iterations = NULL, remove_border_first = TRUE)
}
\arguments{
  \item{data_dir}{Visium data directory.}
  \item{output_dir}{Directory to save tissue_positions snapshots.}
  \item{p_threshold}{Stop when edge vs interior p-value >= threshold.}
  \item{initial_edge_distance}{Starting edge distance.}
  \item{max_iterations}{Optional maximum iterations.}
  \item{remove_border_first}{Whether to remove border before edge testing.}
}
\value{A data.frame with edge distances and p-values.}
\description{Iteratively removes edge layers until the edge matches interior tissue.}
