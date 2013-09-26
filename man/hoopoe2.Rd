\name{hoopoe2}
\alias{hoopoe2}
\title{
Sunrise and sunset times: From light intensity measurement (hoopoe1)
}
\description{
Sunrise and sunset times derived from light intensity measurements over time (\code{\link{hoopoe1}}). The light measurements corresponding to the first part of the annual migration of a European Hoopoe (\emph{Upupa epops}).
}
\usage{data(hoopoe1)}
\format{
A table with 340 rows and 3 columns. Each row corresponds to subsequent twilight events ("tFirst" and "tSecond"). The third column ("type") indicates weather the first event is sunrise (1) or sunset (2). All dates/times are measured in Universal Time Zone (UTC).
}
\source{
Baechler, E., Hahn, S., Schaub, M., Arlettaz, R., Jenni, L., Fox, J.W., Afanasyev, V. & Liechti, F. (2010) Year-Round Tracking of Small Trans-Saharan Migrants Using Light-Level Geolocators. \emph{Plos One}, \bold{5}.
}
\examples{
data(hoopoe2)
attach(hoopoe2)
coord <- coord(tFirst,tSecond,type,degElevation=-6)
## plot in a map using package maps
# par(oma=c(5,0,0,0))
# map(xlim=c(-20,40),ylim=c(-10,60),interior=F,col="darkgrey")
# map(xlim=c(-20,40),ylim=c(-10,60),boundary=F,lty=2,col="darkgrey",add=T)
# mtext(c("Longitude (degrees)","Latitude (degrees)"),side=c(1,2),line=c(2.2,2.5),font=3)
# map.axes()
# points(coord,col="brown",cex=0.5,pch=20)
}
