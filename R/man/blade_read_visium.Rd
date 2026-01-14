\name{blade_read_visium}
\alias{blade_read_visium}
\title{Read standard Visium inputs}
\usage{
blade_read_visium(dir, h5_path = NULL, neighbor_mode = "radius", k_neighbors = 6, neighbor_radius = NULL)
}
\arguments{
  \item{dir}{Path to a Spaceranger output folder.}
  \item{h5_path}{Optional path to filtered_feature_bc_matrix.h5.}
  \item{neighbor_mode}{"radius" or "k" for neighbor selection.}
  \item{k_neighbors}{Number of neighbors for k mode.}
  \item{neighbor_radius}{Optional radius for neighbor search.}
}
\value{A list containing matrix, positions, inclusion mask, and neighbor index.}
\description{Loads Visium matrix/positions and prepares a neighbor index.}
