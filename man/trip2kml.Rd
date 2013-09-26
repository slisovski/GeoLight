\name{trip2kml}
\alias{trip2kml}
\title{Write a file which plots a trip in Google Earth}
\description{This function creates a .kml file from light intensity measurements over time that can ve viewed as a trip in Google Earth.}
\usage{trip2kml(file, tFirst, tSecond, type, degElevation, col.scheme="heat.colors", point.alpha=0.7, cex=1, line.col="goldenrod") }
\arguments{
  \item{file}{A character expression giving the whole path and the name of the resulting output file including the .kml extension.}
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively}
  \item{degElevation}{sun elevation angle in degrees (e.g. -6 for "civil twilight"). Either a single value, a \code{vector} with the same length as \code{tFirst}.}
  \item{col.scheme}{the color scheme used for the points. Possible color schemes are: \code{\link{rainbow}}, \code{\link{heat.colors}}, \code{\link{topo.colors}}, \code{\link{terrain.colors}}.}
  \item{point.alpha}{a \code{numerical value} indicating the transparency of the point colors on a scale from 0 (transparent) to 1 (opaque).}
  \item{cex}{\code{numerical value} for the size of the points.}
  \item{line.col}{An character expression (any of \code{\link{colors}} or hexadecimal notation), or numeric indicating the color of the line connecting the point locations.}
}
\value{This function returns no data. It creates a .kml file in the in the defined path.}
\author{Simeon Lisovski and Michael U. Kemp}
\examples{
data(hoopoe2)
attach(hoopoe2)
filter <- distanceFilter(tFirst,tSecond,type,distance=30)
trip2kml("trip.kml", tFirst[filter], tSecond[filter], type[filter], 
		degElevation=-6, col.scheme="heat.colors", cex=0.7, 
		line.col="goldenrod")
}