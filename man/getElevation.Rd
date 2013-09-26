\name{getElevation}
\alias{getElevation}
\title{Calculate the appropriate sun elevation angle for known location}
\description{Function to calculate the sun elevation angle for light measurements at a known location}
\usage{getElevation(tFirst,tSecond,type,known.coord,plot=TRUE)
}

\arguments{
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining tFirst as sunrise or sunset respectively}
  \item{known.coord}{a \code{SpatialPoint} or \code{matrix} object, containing known x and y coordinates (in that order) for the selected measurement period.}
  \item{plot}{\code{logical}, if TRUE a plot will be produced.}
}

\references{
Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and precision affected by environmental factors. \emph{Methods in Ecology and Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
}
\author{Simeon Lisovski}
\examples{
data(calib2)
attach(calib2)
known.coord <- c(7.1,46.3)
getElevation(tFirst,tSecond,type,known.coord)
}