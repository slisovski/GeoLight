% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/GeoLight.R
\name{getElevation}
\alias{getElevation}
\title{Calculate the appropriate sun elevation angle for known location}
\usage{
getElevation(tFirst, tSecond, type, twl, known.coord, plot = TRUE,
  lnorm.pars = FALSE)
}
\arguments{
\item{tFirst}{vector of sunrise/sunset times (e.g. 2008-12-01 08:30).}

\item{tSecond}{vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).}

\item{type}{vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.}

\item{twl}{data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).}

\item{known.coord}{a \code{SpatialPoint} or \code{matrix} object, containing
known x and y coordinates (in that order) for the selected measurement
period.}

\item{plot}{\code{logical}, if TRUE a plot will be produced.}

\item{lnorm.pars}{\code{logical}, if TRUE shape and scale parameters of the twilight error (log-normal distribution)
will be estimated and included in the output (see Details).}
}
\description{
Function to calculate the median sun elevation angle for light measurements at a
known location and the choosen light threshold.
}
\details{
Optionally, shape and scale paramters of the twiligth error (in minutes) can be estimated. The error is assumed
to follow a log-normal distribution and 0 (elev0) is set 0.1 below the minimum sun elevation angle of estimated twilight times.
Those parameters might be of interest for sensitivity analysis or further processing using the R Package SGAT (https://github.com/SWotherspoon/SGAT).
}
\examples{
data(calib2)
  calib2$tFirst  <- as.POSIXct(calib2$tFirst, tz = "GMT")
  calib2$tSecond <- as.POSIXct(calib2$tSecond, tz = "GMT")
getElevation(calib2, known.coord = c(7.1,46.3), lnorm.pars = TRUE)
}
\author{
Simeon Lisovski
}
\references{
Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt,
F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
precision affected by environmental factors. \emph{Methods in Ecology and
Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
}

