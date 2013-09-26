\name{coord}
\alias{coord}
\title{Threshold based geographical positioning}
\description{Calculate the latitude and longitude from two subsequent twilight events (sunrise/sunset)}
\usage{coord(tFirst,tSecond,type,degElevation=-6,site=FALSE,sites,note=TRUE)}
\arguments{
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively}
  \item{degElevation}{sun elevation angle in degrees (e.g. -6 for "civil twilight"). Either a single value, a \code{vector} with the same length as \code{tFirst}. If site=TRUE a \code{vector} with a sun elevation angle for each site in numerical order of the sites (incl. 0 for movements).}
  \item{site}{\code{logical},if TRUE the a sun elevation angle for each site (incl. 0 dor movements) will be used for positioning}
  \item{sites}{a \code{numerical vector} assigning each row to a particular period. Stationary periods in numerical order and values >0, migration/movement periods 0.}
	\item{note}{\code{logical}, if TRUE a notation of how many positions could be calculated in proportion to the number of failures will be printed at the end.}
}

\details{The format (date and time) of \emph{tFirst} and \emph{tSecond} has to be "yyyy-mm-dd hh:mm" corresponding to Universal Time Zone UTC (see: \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones}).}
\value{A matrix of coordinates in decimal degrees. First column are longitudes, expressed in degrees east of Greenwich. Second column contains the latitudes in degrees north the equator. If latitude can not be calculated (e.g. during equinox, at high latitudes) the function will return NA for the date.
}
\references{Montenbruck, O. & Pfleger, T. (2000) Astronomy on the Personal Computer. \emph{Springer}, Berlin.}
\author{Simeon Lisovski}
\examples{
data(hoopoe2)
attach(hoopoe2)
coord <- coord(tFirst,tSecond,type,degElevation=-6)
## tripMap(coord,xlim=c(-20,20),ylim=c(5,50), main="hoopoe2")
}