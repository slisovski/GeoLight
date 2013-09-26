
#' The GeoLight Package
#'
#' This is a summary of all features of \bold{\code{GeoLight}}, a
#' \code{R}-package for analyzing light based geolocator data.
#'
#' @name GeoLIght-package
#' @docType package
#' @author Simeon Lisovski, Silke Bauer, Tamara Emmenegger
#' @aliases GeoLight
#' @section Details: \bold{\code{GeoLight}} is a package to derive geographical
#' positions from daily light intensity pattern. Positioning and calibration
#' methods are based on the threshold-method (Ekstrom 2004, Lisovski \emph{et
#' al.} 2012). A changepoint model from the \code{R} package \code{changepoint}
#' is implemented to distinguish between periods of residency and movement
#' based on the sunrise and sunset times. Mapping functions are implemented
#' using the \code{R} package \code{maps}.
NULL



#' Example data for calibration: Light intensities and twilight events
#'
#' Light intensity measurements over time (calib1) recorded at the rooftop of
#' the Swiss Ornithological Institute (Lon: 8.0, Lat: 47.01). Defined twilight
#' events from calib1 (calib2). These data serve as an example for calculating
#' the sun elevation angle of an additional data set, which is subsequently
#' used to calibrate the focal dataset.
#'
#' @name calib1
#' @docType data
#' @aliases calib1 calib2 calib
#' @references Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt,
#' F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
#' precision affected by environmental factors. \emph{Methods in Ecology and
#' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
#' @examples
#'
#' data(calib2)
#' attach(calib2)
#' getElevation(tFirst,tSecond,type,c(8,47.01))
#'
NULL


#' Light intensity measurements over time recorded on a migratory bird
#'
#' Sunlight intensity measurements over time recorded during the first part of
#' the annual migration of a European Hoopoe (\cite{Upupa epops}). All
#' dates/times are measured in Universal Time Zone (UTC).
#'
#'
#' @name hoopoe1
#' @docType data
#' @format A table with 24474 rows and 2 columns, rows corresponding to light
#' measurements recorded in ten-minute intervals (datetime).
#' @source Baechler, E., Hahn, S., Schaub, M., Arlettaz, R., Jenni, L., Fox,
#' J.W., Afanasyev, V. & Liechti, F. (2010) Year-Round Tracking of Small
#' Trans-Saharan Migrants Using Light-Level Geolocators. \emph{Plos One},
#' \bold{5}.
NULL


#' Sunrise and sunset times: From light intensity measurement (hoopoe1)
#'
#' Sunrise and sunset times derived from light intensity measurements over time
#' (\code{\link{hoopoe1}}). The light measurements corresponding to the first
#' part of the annual migration of a European Hoopoe (\emph{Upupa epops}).
#'
#' @name hoopoe2
#' @docType data
#' @format A table with 340 rows and 3 columns. Each row corresponds to
#' subsequent twilight events ("tFirst" and "tSecond"). The third column
#' ("type") indicates weather the first event is sunrise (1) or sunset (2). All
#' dates/times are measured in Universal Time Zone (UTC).
#' @source Baechler, E., Hahn, S., Schaub, M., Arlettaz, R., Jenni, L., Fox,
#' J.W., Afanasyev, V. & Liechti, F. (2010) Year-Round Tracking of Small
#' Trans-Saharan Migrants Using Light-Level Geolocators. \emph{Plos One},
#' \bold{5}.
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' coord <- coord(tFirst,tSecond,type,degElevation=-6)
#' ## plot in a map using package maps
#' # par(oma=c(5,0,0,0))
#' # map(xlim=c(-20,40),ylim=c(-10,60),interior=F,col="darkgrey")
#' # map(xlim=c(-20,40),ylim=c(-10,60),boundary=F,lty=2,col="darkgrey",add=T)
#' # mtext(c("Longitude (degrees)","Latitude (degrees)"),side=c(1,2),line=c(2.2,2.5),font=3)
#' # map.axes()
#' # points(coord,col="brown",cex=0.5,pch=20)
#'
NULL



