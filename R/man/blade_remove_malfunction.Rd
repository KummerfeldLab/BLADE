\name{blade_remove_malfunction}
\alias{blade_remove_malfunction}
\title{Remove malfunction (outlier) spots}
\usage{
blade_remove_malfunction(tissue, out_rate = 0.05, use_clusters = TRUE,
  min_cluster_size = 15, radius = NULL, radius_scale = 1.0)
}
\arguments{
  \item{tissue}{A BLADE tissue list.}
  \item{out_rate}{Expected outlier rate.}
  \item{use_clusters}{If TRUE, only removes large outlier clusters.}
  \item{min_cluster_size}{Minimum cluster size to remove.}
  \item{radius}{Neighbor radius for clustering.}
  \item{radius_scale}{Scale factor for the radius.}
}
\value{Updated tissue list.}
\description{Marks outlier spots as excluded, optionally clustering first.}
