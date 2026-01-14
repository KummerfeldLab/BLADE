\name{blade_remove_border}
\alias{blade_remove_border}
\title{Remove border spots}
\usage{
blade_remove_border(tissue, distance = 0.0, strategy = "extreme", expected_neighbors = 8)
}
\arguments{
  \item{tissue}{A BLADE tissue list.}
  \item{distance}{Border width for extreme strategy.}
  \item{strategy}{"extreme" or "degree".}
  \item{expected_neighbors}{Neighbor threshold for degree strategy.}
}
\value{Updated tissue list.}
\description{Marks border spots as excluded.}
