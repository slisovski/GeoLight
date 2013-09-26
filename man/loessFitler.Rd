\name{loessFilter}
\alias{loessFilter}
\title{Filter to remove outliers in defined sunrise and sunset times}
\description{This filter defines outliers based on residuals from a local polynomial regression fitting provcess (\code{\link{loess}}).}
\usage{loessFilter(tFirst, tSecond, type, k=3, plot=TRUE)}
\arguments{
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively}
  \item{k}{a measure of how many interquartile ranges to take before saying that a particular twilight event is an outlier}
  \item{plot}{code{logical}, if TRUE a plot indicating the filtered times will be produced.}
}
\value{Logical \code{vector} matching positions that pass the filter.}
\author{Simeon Lisovski & Eldar Rakhimberdiev}