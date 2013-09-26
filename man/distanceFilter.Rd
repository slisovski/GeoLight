\name{distanceFilter}
\alias{distanceFilter}
\title{Filter for unrealistic positions within a track}
\description{The filter identifies unrealistic positions. The maximal distance per hour/day can be set corresponding to the particular species.}
\usage{distanceFilter(tFirst,tSecond,type,degElevation,distance,units)}
\arguments{
   \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining tFirst as sunrise or sunset respectively}
  \item{degElevation}{sun elevation angle in degrees (e.g. -6 for "civil twilight")}
  \item{distance}{the maximal distance in km per \code{units}. Distances above will be considered as unrealistic.}
  \item{units}{the time unite corresponding to the distance. Default is "hour", alternative option is "day".}
}
\value{Logical \code{vector} matching positions that pass the filter.}
\details{The filter currently depends to some extent on the calibration (\code{degElevation}).}
\author{Simeon Lisovski, Fraenzi Korner-Nievergelt}
\examples{
data(hoopoe2)
attach(hoopoe2)
coord  <- coord(tFirst,tSecond,type)
filter <- distanceFilter(tFirst,tSecond,type,distance=30)
tripMap(coord[filter,],xlim=c(-20,20),ylim=c(0,60),main="hoopoe2 (filter)")
}