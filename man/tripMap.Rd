\name{tripMap}
\alias{tripMap}
\title{Draw the positions and the trip on a map}
\description{Draw a map (from the \code{R} Package \code{maps}) with calculated positions connected by a line}
\usage{tripMap(coord,equinox,
	map.range=c("EuroAfrica","AustralAsia","America","World"),
	xlim,ylim,xlab,ylab,col,main, legend,add)}
\arguments{
  \item{coord}{a \code{SpatialPoints} or \code{matrix} object, containing x and y coordinates (in that order).}
  \item{equinox}{logical; if \code{TRUE}, the equinox period(s) is shown as a broken blue line.}
  \item{map.range}{some possibilities to choose defined areas (default: "World").}
\item{xlim}{two element numeric vector giving a range of longitudes, expressed in degrees, to which the map is restricted. Longitude is measured in degrees east of Greenwich, so that, in particular, locations in Switzerland have positive longitude. If map.range is defined this argument will not be considered.}
  \item{ylim}{two element numeric vector giving a range of latitudes, expressed in degrees, to which the map is restricted. Latitude is measured in degrees north of the equator, so that, in particular, locations in Switzerland have positive latitude. If map.range is defined this argument will not be considered.}
 \item{xlab}{a title for the x axis.} 
 \item{ylab}{a title for the y axis.} 
 \item{col}{colour code or name for the trip line (default equals red).} 
 \item{main}{map title}
 \item{legend}{\code{logical}; if \code{TRUE}, a legend will be added to the plot.}
 \item{add}{\code{logical}; if \code{TRUE}, positions will be added to an existing plot.}
}
\author{Simeon Lisovski}
\examples{
data(hoopoe2)
attach(hoopoe2)
coord <- coord(tFirst,tSecond,type,-6)
tripMap(coord,xlim=c(-20,20),ylim=c(0,60),main="hoopoe2")
}