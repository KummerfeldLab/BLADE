\name{blade_read_visium_hd}
\alias{blade_read_visium_hd}
\title{Read Visium HD inputs}
\usage{
blade_read_visium_hd(data_dir = NULL, h5_path = NULL, positions_path = NULL,
  neighbor_mode = "radius", k_neighbors = 8, neighbor_radius = NULL)
}
\arguments{
  \item{data_dir}{Optional data directory for auto-discovery.}
  \item{h5_path}{Path to a 10x HDF5 matrix file.}
  \item{positions_path}{Path to tissue_positions.parquet or CSV.}
  \item{neighbor_mode}{"radius" or "k" for neighbor selection.}
  \item{k_neighbors}{Number of neighbors for k mode.}
  \item{neighbor_radius}{Optional radius for neighbor search.}
}
\value{A list containing matrix, positions, inclusion mask, and neighbor index.}
\description{Loads Visium HD HDF5 inputs and positions with a neighbor index.}
