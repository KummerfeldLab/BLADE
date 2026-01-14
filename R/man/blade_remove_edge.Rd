\name{blade_remove_edge}
\alias{blade_remove_edge}
\title{Remove edge spots}
\usage{
blade_remove_edge(tissue, distance = 1, strategy = "extreme", expected_neighbors = 8,
  seed_mode = "zeros", border_width = NULL)
}
\arguments{
  \item{tissue}{A BLADE tissue list.}
  \item{distance}{Edge depth to remove.}
  \item{strategy}{"extreme" or "degree" for border seeding.}
  \item{expected_neighbors}{Neighbor threshold for degree strategy.}
  \item{seed_mode}{"zeros", "border", or "both".}
  \item{border_width}{Border width for extreme strategy.}
}
\value{Updated tissue list.}
\description{Marks edge spots up to the requested depth as excluded.}
