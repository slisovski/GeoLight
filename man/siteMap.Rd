\name{siteMap}
\alias{siteMap}
\title{Draws sites of residency and adds a convex hull}
\description{Draw a map (from the \code{R} Package \code{maps}) showing the defined stationary sites}
\usage{siteMap(coord,site,points=TRUE,
	map.range=c("EuroAfrica","AustralAsia","America","World"),
	xlim,ylim,xlab="Longitude",ylab="Latitude",lwd=1,lty=1,pch=1,cex=1,col,main,add=FALSE)}
\arguments{
  \item{coord}{a \code{SpatialPoints} or \code{matrix} object, containing x and y coordinates (in that order).}
  \item{site}{a \code{numerical vector} assigning each row to a particular period. Stationary periods in numerical order and values >0, migration/movement periods 0.}
  \item{points}{\code{logical}; if \code{TRUE}, the points of each site will also be plottet.}
  \item{map.range}{some possibilities to choose defined areas ("World (default)", "EuroAfrica","America","AustralAsia").}
  \item{xlim}{two element numeric vector giving a range of longitudes, expressed in degrees, to which the map is restricted. Longitude is measured in degrees east of Greenwich, so that, in particular, locations in Switzerland have positive longitude. If map.range is defined this argument will not be considered.}
  \item{ylim}{two element numeric vector giving a range of latitudes, expressed in degrees, to which the map is restricted. Latitude is measured in degrees north of the equator, so that, in particular, locations in Switzerland have positive latitude. If map.range is defined this argument will not be considered.}
 \item{xlab}{a title for the x axis.} 
 \item{ylab}{a title for the y axis.}
 \item{lwd}{The line width, a positive number, defaulting to 1.} 
 \item{lty}{The line type. Line types can either be specified as an integer (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) or as one of the character strings "blank", "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash", where "blank" uses "invisible lines" (i.e., does not draw them).}
 \item{pch}{Either an integer specifying a symbol or a single character to be used as the default in plotting points. See \code{\link{points}} for possible values and their interpretation. Note that only integers and single-character strings can be set as a graphics parameter (and not NA nor NULL).}   
 \item{cex}{A numerical value giving the amount by which plotting symbols should be magnified relative to the default.} 
 \item{col}{a vector of colors with the same length as the number of defined sites (if Default a predefined color ramp will be used).} 
 \item{main}{map title.}
 \item{add}{\code{logical}; if \code{TRUE}, positions will be added to an existing plot.}
}
\author{Simeon Lisovski}
\examples{
data(hoopoe2)
attach(hoopoe2)
coord <- coord(tFirst,tSecond,type,-6)
filter <- distanceFilter(tFirst,tSecond,type,distance=30)
site <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)$site
siteMap(coord[filter,],site[filter],xlim=c(-20,20),ylim=c(0,60),lwd=2,pch=20,cex=0.5,main="hoopoe2")
}