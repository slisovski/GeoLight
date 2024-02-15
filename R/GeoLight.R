#' The GeoLight Package
#' @aliases GeoLight
#' @aliases GeoLight-package
#' @description This is a summary of all features of \bold{\code{GeoLight}}, a \code{R}-Package For 
#' Analyzing Light Based Geolocator Data
#' @details \bold{\code{GeoLight}} is a package to derive geographical positions from daily light intensity pattern. 
#' Positioning and calibration methods are based on the threshold-method (Ekstrom 2004, Lisovski \emph{et al.} 2012). 
#' A changepoint model from the \code{R} package \code{changepoint} is implemented to distinguish between periods of 
#' residency and movement based on the sunrise and sunset times. Mapping functions are implemented 
#' using the \code{R} package \code{maps}.
#' @section Getting Started: 
#' We refrain from giving detailed background on the (several steps of) 
#' analysis of light-based Geolocator data here but strongly recommend the key-publications below. 
#' @section Updates:
#' We advise all users to update their installation of \bold{\code{GeoLight}} regularly.
#' Type \code{news(package="GeoLight")} to read news documentation about changes to the recent and all previous version of the package
#' @section Important notes:
#' Most functions in \bold{\code{GeoLight}} require the same initial units and mostly the format and object type is mandatory:
#' \tabular{rll}{
#' \tab \code{tFirst} \tab yyyy-mm-dd hh:mm "UTC" (see: \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones})\cr
#' \tab \code{tSecond} \tab as \emph{tFirst} (e.g. 2008-12-01 17:30) \cr
#' \tab \code{type} \tab either 1 or 2 depending on wheter \emph{tFirst} is sunrise (1) or sunset (2)\cr
#' \tab \code{coord} \tab \code{SpatialPoints} or a \code{matrix}, containing x and y coordinates (in that order) \cr
#' \tab \code{degElevation} \tab a \code{vector} or a single \code{value} of sun elevation angle(s) in degrees (e.g. -6)
#' }
#' @section FUNCTIONS AND DATASETS:
#' In the following, we give a summary of the main functions and sample datasets in the \bold{\code{GeoLight}} package. 
#' Alternatively a list of all functions and datasets in alphabetical order is available by typing \code{library(help=GeoLight)}. 
#' For further information on any of these functions, type \code{help(function name)}.
#'
#'  @section CONTENTS:
#' \tabular{rll}{
#' \tab {I.} \tab Determination of sunset and sunrise \cr
#' \tab {II.} \tab Residency analysis \cr
#' \tab {III.} \tab Calibration \cr
#' \tab {IV.} \tab {Positioning}\cr
#' \tab {V.} \tab {Data visualisation}\cr
#' \tab {VI.} \tab {Examples}
#' }
#' @section I. Determination of sunset and sunrise:
#' \tabular{rll}{
#' \tab \code{\link{gleTrans}} \tab transformation of already defined twilight events* \cr
#' \tab \code{\link{glfTrans}} \tab transformation of light intensity measurements over time* \cr
#' \tab \code{\link{luxTrans}} \tab transformation of light intensity measurements over time** \cr
#' \tab \code{\link{lightFilter}} \tab filter to remove noise in light intensity measurements during the night \cr
#' \tab \code{\link{twilightCalc}} \tab definition of twilight events (\emph{sunrise, sunset}) from light intensity measurements \cr
#' }
#' 
#' * written for data recorded by geolocator devices from the \bold{Swiss Ornithological Institute} \cr
#' ** written for data recorded by geolocator devices from \bold{Migrate Technology Ltd}
#' 
#' @section II. Residency Analysis:
#' \tabular{rll}{
#' \tab \code{\link{changeLight}} \tab function to distinguish between residency and movement periods \cr
#' \tab \code{\link{schedule}} \tab function to produce a data frame summerizing the residency and movement pattern \cr
#' }
#' 
#' @section III. Calibration:
#' See Lisovski \emph{et al.} 2012 for all implemented calibration methods.
#' \tabular{rll}{
#' \tab \code{\link{getElevation}} \tab function to calculate the sun elevation angle for data with known position \cr
#' \tab \code{\link{HillEkstromCalib}} \tab \emph{Hill-Ekstrom calibration} for one or more defined stationary periods \cr
#' }
#' 
#' @section IV. Positioning:
#' \tabular{rll}{
#' \tab \code{\link{coord}} \tab main function to derive a \code{matrix} of spatial coordinates \cr
#' \tab \code{\link{distanceFilter}} \tab filter function to reduce unrealistic positions (not recommended, since the filtering ignore positioning error) \cr
#' \tab \code{\link{loessFilter}} \tab filter function to define outliers in sunrise and sunset times (defined twilight events) \cr
#' }
#' 
#' @section V. Data visualisation:
#' \tabular{rll}{
#' \tab \code{\link{tripMap}} \tab function to map the derived positions and combine the coordinates in time order\cr
#' \tab \code{\link{siteMap}} \tab function to show the results of the residency analysis on a map \cr
#' }
#' 
#' @section IV. Examples:
#' \tabular{rll}{
#' \tab \code{\link{calib1}} \tab data for calibration: light intensities \cr
#' \tab \code{\link{calib2}} \tab  data for calibration: Calculated twilight events (from \code{\link{calib1}} by \code{\link{twilightCalc}}) \cr
#' \tab \code{\link{hoopoe1}} \tab light intensity measurements over time recorded on a migratory bird \cr
#' \tab \code{\link{hoopoe2}} \tab sunrise and sunset times: From light intensity measurement (from \code{\link{hoopoe1}}) \cr
#' }
#' 
#' @section \bold{\code{R}} Packages for Further Spatial Ananlyses:
#' \code{spatstat} \cr
#' \code{adehabitat} \cr
#' \code{gstat} \cr
#' \code{trip} \cr
#' \code{tripEstimation} \cr
#' \code{move} \cr
#' ...
#' 
#' @section Acknowledgements:
#' Steffen Hahn, Felix Liechti, Fraenzi Korner-Nievergelt, Andrea Koelzsch, Eldar Rakhimberdiev, Erich Baechler, Eli Bridge,  Andrew Parnell, Richard Inger
#'
#' @section Authors:
#' Simeon Lisovski, Simon Wotherspoon, Michael Sumner, Silke Bauer, Tamara Emmenegger \cr
#' Maintainer: Simeon Lisovski <simeon.lisovski(at)gmail.com>
#' 
#' @section References:
#' Ekstrom, P.A. (2004) An advance in geolocation by light. \emph{Memoirs of the National Institute of Polar Research}, Special Issue, \bold{58}, 210-226.
#' 
#' Fudickar, A.M., Wikelski, M., Partecke, J. (2011) Tracking migratory songbirds: accuracy of light-level loggers (geolocators) in forest habitats. \emph{Methods in Ecology and Evolution}, DOI: 10.1111/j.2041-210X.2011.00136.x.
#'
#' Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step: Latitude. \emph{Electronic Tagging and Tracking in Marine Fisheries} (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The Netherlands.
#'
#' Hill, R.D. (1994) Theory of geolocation by light levels. \emph{Elephant Seals: Population Ecology, Behavior, and Physiology} (eds L. Boeuf, J. Burney & R.M. Laws), pp. 228-237. University of California Press, Berkeley.
#'  
#' Lisovski, S. and Hahn, S. (2012) GeoLight - processing and analysing light-based geolocator data in R. \emph{Methods in Ecology and Evolution}, doi: 10.1111/j.2041-210X.2012.00248.x
#' 
#' Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and precision affected by environmental factors. \emph{Methods in Ecology and Evolution}, doi: 10.1111/j.2041-210X.2012.00185.x
#' 
#' Wilson, R.P., Ducamp, J.J., Rees, G., Culik, B.M. & Niekamp, K. (1992) Estimation of location: global coverage using light intensity. \emph{Wildlife telemetry: remote monitoring and tracking of animals} (eds I.M. Priede & S.M. Swift), pp. 131-134. Ellis Horward, Chichester.

##' @name globalVariables
##' @import utils
globalVariables(".")

##' @title Checkup of arguments in GeoLight tables
##' @usage i.argCheck(y)
##' @param y GeoLigth data table.
##' @author Simeon Lisovski
##' @noRd
i.argCheck <- function(y) {
  if(any(sapply(y, function(x) class(x))=="data.frame")) {
    ind01 <- which(sapply(y, function(x) class(x))=="data.frame")
    if(!all(ind02 <- c("tFirst", "tSecond", "type")%in%names(y[[ind01]]))) {
      whc <- paste("The following columns in data frame twl are missing with no default: ", paste(c("tFirst", "tSecond", "type")[!ind02], collapse = " and "), sep = "")
      stop(whc , call. = F)
    } 
    out <- y[[ind01]]
  } else {
    if(!all(c("tFirst", "tSecond", "type")%in%names(y))) {
      ind03 <- c("tFirst", "tSecond", "type")%in%names(y)
      stop(sprintf(paste(paste(c("tFirst", "tSecond", "type")[!ind03], collapse = " and "), "is missing with no default.")))
    } else {
      out <- data.frame(tFirst = y$tFirst, tSecond = y$tSecond, type = y$type)
    }
  }
  if(any(c(class(out[,1])[1], class(out[,2])[1])!="POSIXct")) {
    stop(sprintf("Date and time inforamtion (e.g. tFirst and tSecond) need to be provided as POSIXct class objects."), call. = F)
  }
  out  
}

##' Estimate location from consecutive twilights

##' This function estimates the location given the times at which 
##' the observer sees two successive twilights.
##' 
##' Longitude is estimated by computing apparent time of local noon
##' from sunrise and sunset, and determining the longitude for which
##' this is noon. Latitude is estimated from the required zenith and
##' the sun's hour angle for both sunrise and sunset, and averaged.
##'
##' When the solar declination is near zero (at the equinoxes)
##' latitude estimates are extremely sensitive to errors.  Where the
##' sine of the solar declination is less than \code{tol}, the
##' latitude estimates are returned as \code{NA}.
##' 
##' The format (date and time) of \emph{tFirst} and \emph{tSecond} has to be
##' "yyyy-mm-dd hh:mm" corresponding to Universal Time Zone UTC (see:
##' \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones})
##' 
##' @title Simple Threshold Geolocation Estimates
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param degElevation the sun elevation angle (in degrees) that defines twilight (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst} or \code{nrow(x)}.
##' @param tol tolerance on the sine of the solar declination (only implemented in method 'NOAA').
##' @param note \code{logical}, if TRUE a notation of how many positions could
##' be calculated in proportion to the number of failures will be printed at the
##' end.
##' @param method Defines the method for the location estimates. 'NOAA' is based on
##' code and the excel spreadsheet from the NOAA site (http://www.esrl.noaa.gov/gmd/grad/solcalc/),
##' 'Montenbruck' is based on Montenbruck, O. & Pfleger, T. (2000) Astronomy on the Personal
##' Computer. \emph{Springer}, Berlin.
##' @return A matrix of coordinates in decimal degrees. First column are
##' longitudes, expressed in degrees east of Greenwich. Second column contains
##' the latitudes in degrees north the equator.
##' @author Simeon Lisovski, Simon Wotherspoon, Michael Sumner
##' @examples
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' crds <- coord(hoopoe2, degElevation=-6, tol = 0.2)
##' tripMap(crds, xlim=c(-20,20), ylim=c(5,50), main="hoopoe2")
##' @export   
coord  <- function(tFirst, tSecond, type, twl, degElevation = -6, tol = 0, method = "NOAA",  note = TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  rise <- ifelse(tab$type==1, tab$tFirst, tab$tSecond)
  set <- ifelse(tab$type==1, tab$tSecond, tab$tFirst)
  
  if(method == "NOAA") {
    rad <- pi/180
    sr <- solar(rise)
    ss <- solar(set)
    cosz <- cos(rad*(90-degElevation))
    lon <- -(sr$solarTime+ss$solarTime+ifelse(sr$solarTime<ss$solarTime,360,0))/2
    lon <- (lon+180)%%360-180
    
    ## Compute latitude from sunrise
    hourAngle <- sr$solarTime+lon-180
    a <- sr$sinSolarDec
    b <- sr$cosSolarDec*cos(rad*hourAngle)
    x <- (a*cosz-sign(a)*b*suppressWarnings(sqrt(a^2+b^2-cosz^2)))/(a^2+b^2)
    lat1 <- ifelse(abs(a)>tol,asin(x)/rad,NA)
    
    ## Compute latitude from sunset
    hourAngle <- ss$solarTime+lon-180
    a <- ss$sinSolarDec
    b <- ss$cosSolarDec*cos(rad*hourAngle)
    x <- (a*cosz-sign(a)*b*suppressWarnings(sqrt(a^2+b^2-cosz^2)))/(a^2+b^2)
    lat2 <- ifelse(abs(a)>tol,asin(x)/rad,NA)
    
    ## Average latitudes
    out <- cbind(lon=lon,lat=rowMeans(cbind(lat1,lat2),na.rm=TRUE))
  } 
  
  if(method == "Montenbruck") {
    
    out <- coord2(tab$tFirst, tab$tSecond, tab$type, degElevation)
    
  }
  
  
  if(note) cat(paste("Note: Out of ", nrow(out)," twilight pairs, the calculation of ", sum(is.na(out[,2]))," latitudes failed ","(",
                     floor(sum(is.na(out[,2])*100)/nrow(out))," %)",sep=""))
  out
}

coord2 <- function(tFirst, tSecond, type, degElevation=-6) {
  
  # if noon, RadHourAngle = 0, if midnight RadHourAngle = pi
  # --------------------------------------------------------
  RadHourAngle <- numeric(length(type))
  index1 <- type==1
  if (sum(index1)>0) RadHourAngle[index1] <- 0
  RadHourAngle[!index1] <- pi
  # --------------------------------------------------------
  
  tSunTransit <- as.character(as.POSIXct(tFirst, tz = "GMT") + as.numeric(difftime(as.POSIXct(tSecond, tz = "GMT"), 
                                                                                   as.POSIXct(tFirst, tz = "GMT"), 
                                                                                   units="secs")/2))
  
  index0 <- (nchar(tSunTransit) <= 10)
  if(sum(index0)>0) tSunTransit[index0] <- as.character(paste(tSunTransit[index0]," ","00:00",sep=""))
  
  # Longitude
  jD  <- i.julianDate(as.numeric(substring(tSunTransit,1,4)),as.numeric(substring(tSunTransit,6,7)),
                      as.numeric(substring(tSunTransit,9,10)),as.numeric(substring(tSunTransit,12,13)),as.numeric(substring(tSunTransit,15,16)))
  jD0 <- i.julianDate(as.numeric(substring(tSunTransit,1,4)),as.numeric(substring(tSunTransit,6,7)),
                      as.numeric(substring(tSunTransit,9,10)),rep(0,length(tFirst)),rep(0,length(tFirst)))
  jC  <- i.JC2000(jD)
  jC0 <- i.JC2000(jD0)
  
  radObliquity         <- i.radObliquity(jC)
  radEclipticLongitude <- i.radEclipticLongitude(jC)
  radRightAscension    <- i.radRightAscension(radEclipticLongitude,radObliquity)
  radGMST              <- i.radGMST(jD,jD0,jC,jC0)
  
  degLongitude <- i.setToRange(-180,180,i.deg(RadHourAngle + radRightAscension - radGMST))
  
  
  # Latitude
  jD  <- i.julianDate(as.numeric(substring(tFirst,1,4)),as.numeric(substring(tFirst,6,7)),
                      as.numeric(substring(tFirst,9,10)),as.numeric(substring(tFirst,12,13)),
                      as.numeric(substring(tFirst,15,16)))
  jD0 <- i.julianDate(as.numeric(substring(tFirst,1,4)),as.numeric(substring(tFirst,6,7)),
                      as.numeric(substring(tFirst,9,10)),rep(0,length(tFirst)),rep(0,length(tSecond)))
  jC  <- i.JC2000(jD)
  jC0 <- i.JC2000(jD0)
  
  
  radElevation         <- if(length(degElevation)==1) rep(i.rad(degElevation),length(jD)) else i.rad(degElevation)
  
  sinElevation         <- sin(radElevation)
  radObliquity         <- i.radObliquity(jC)
  radEclipticLongitude <- i.radEclipticLongitude(jC)
  radDeclination       <- i.radDeclination(radEclipticLongitude,radObliquity)
  sinDeclination       <- sin(radDeclination)
  cosDeclination       <- cos(radDeclination)
  
  radHourAngle         <- i.radGMST(jD,jD0,jC,jC0) + i.rad(degLongitude) - i.radRightAscension(radEclipticLongitude,radObliquity)
  cosHourAngle         <- cos(radHourAngle)
  
  term1                <- sinElevation/(sqrt(sinDeclination^2 + (cosDeclination*cosHourAngle)^2))
  term2                <- (cosDeclination*cosHourAngle)/sinDeclination
  
  degLatitude <- numeric(length(radElevation))
  
  index1  <- (abs(radElevation) > abs(radDeclination))
  if (sum(index1)>0) degLatitude[index1] <- NA
  
  index2  <- (radDeclination > 0 & !index1)
  if (sum(index2)>0)  degLatitude[index2] <- i.setToRange(-90,90,i.deg(asin(term1[index2])-atan(term2[index2])))
  
  index3  <- (radDeclination < 0 & !index1)
  if (sum(index3)>0)  degLatitude[index3] <- i.setToRange(-90,90,i.deg(acos(term1[index3])+atan(1/term2[index3])))
  
  degLatitude[!index1&!index2&!index3] <- NA
  
  cbind(degLongitude, degLatitude)
}

##' Function to calculate the sun elevation angle for light measurements at a
##' known location and the choosen light threshold.
##' 
##' NEW: The function provides two different sun elevation angle. The first (a1) refers to the median twiligth error and
##' is needed for threshold based locaiton estimation (e.g. \code{\link{coord}}). The second (a0) refers to the zero elevation
##' angle of the twilight error distribution and is required for e.g. \code{\link{mergeSites2}}, \code{\link{siteEstimate}}. The function
##' also provides the parameters (log.mean and log.sd) of the fitted log-normal distribution (see red line in plot),
##'
##' @title Calculate the appropriate sun elevation angle for known location
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param known.coord a \code{SpatialPoint} or \code{matrix} object, containing
##' known x and y coordinates (in that order) for the selected measurement
##' period.
##' @param method the function can either estimate the sun elevation angle and the twilight error parameters using a log-normal ("log-norm")
##' or a gamma ("gamma") error distribution. It is recommended to try both and evaluate the fit using the plot.
##' @param plot \code{logical}, if TRUE a plot will be produced.
##' @author Simeon Lisovski
##' @references Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt,
##' F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
##' precision affected by environmental factors. \emph{Methods in Ecology and
##' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
##' @examples
##' data(calib2)
##'   calib2$tFirst  <- as.POSIXct(calib2$tFirst, tz = "GMT")
##'   calib2$tSecond <- as.POSIXct(calib2$tSecond, tz = "GMT")
##' getElevation(calib2, known.coord = c(7.1,46.3))
##' @importFrom MASS fitdistr
##' @importFrom graphics hist
##' @importFrom stats dlnorm median lm dgamma
##' @importFrom ggplot2 geom_histogram labs theme_minimal geom_point geom_text geom_label annotate
##' @export getElevation

getElevation <- function(tFirst, tSecond, type, twl, known.coord, method = "log-norm", plot=TRUE) {
  
  if(!(method%in%c("gamma", "log-norm"))) stop("Method can only be `gamma` or `log-norm`.")
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   #gets GeoLgiht data
  tab <- geolight.convert(tab[,1], tab[,2], tab[,3])  #converst data to other format where first row beginning time and second row rise or fall
  
  sun  <- solar(tab[,1])  #Calculate solar time, the equation of time and solar declination
  z    <- refracted(zenith(sun, known.coord[1], known.coord[2])) # Adjust the solar zenith angle for atmospheric refraction.
  # plot(z)
  
  inc = 0
  repeat {
    twl_t   <- twilight(tab[,1], known.coord[1], known.coord[2], rise = tab[,2], zenith = max(z)+inc)
    twl_dev <- ifelse(tab$Rise, as.numeric(difftime(tab[,1], twl_t, units = "mins")),
                      as.numeric(difftime(twl_t, tab[,1], units = "mins")))
    if(all(twl_dev>=0)) {
      break
    } else {
      inc <- inc+0.01
    }
  }
  z0 <- max(z)+inc
  
  seq <- seq(0, max(twl_dev), length = 100)
  
  if(method=="log-norm"){
    fitml_ng <- suppressWarnings(fitdistr(twl_dev, "log-normal"))
    lns      <- dlnorm(seq, fitml_ng$estimate[1], fitml_ng$estimate[2])
  }
  if(method=="gamma"){
    fitml_ng <- suppressWarnings(fitdistr(twl_dev, "gamma"))
    lns      <- dgamma(seq, fitml_ng$estimate[1], fitml_ng$estimate[2])
  }
  
  
  diffz <- as.data.frame(cbind(min = apply(cbind(tab[,1], twilight(tab[,1], known.coord[1], known.coord[2], rise = tab[,2], zenith = z0)), 1, function(x) abs(x[1]-x[2]))/60, z = z))
  mod  <- lm(z~min, data = diffz)
  mod2 <- lm(min~z, data = diffz)
  
  a1.0 <- seq[which.max(lns)]
  a1.1 <- 90-predict(mod, newdata = data.frame(min = a1.0))
  
  if(plot) {
    # base histograms
    base_hist <- ggplot() +
      geom_histogram(aes(x = twl_dev),
                     bins = round(max(twl_dev)) + 1 ,
                     #bins =  28,
                     col = "black",
                     fill="lightgrey" )
    
    # Only label changes to log-norm/gamme
    if(method=="log-norm") lab_hist <- base_hist +
        labs(title =  "Twilight Model (log-norm)",
             x = "twilight error (min)",
             y = "Density")
    
    
    if(method=="gamma") lab_hist <- base_hist + 
        labs(title =  "Twilight Model (gamma)",
             x = "twilight error (min)",
             y = "Density")
    
    # base plot + red lines, two points, legends 
    labels <- round(90-predict(mod, newdata = data.frame(min = seq(0, max(twl_dev), 6))),1) # labels for sun elevation degree 
    at <- seq(0, max(twl_dev), 6) # x values for sun elevation degree 
    
    # here we get the warning: "argument ‘freq’ is not made use of", because we use freq and turn of plotting, as freq also adjusts something in the graphic representation. Therefore I silenced the warnings in the following chunk
    height <- rep((max(hist(twl_dev, 
                            breaks = round(max(twl_dev)),
                            plot = FALSE)$counts)) + 1.5,
                  times = 6) # y for sun elevation degree
    
    sun.elv.angle.df <- as.data.frame(cbind(labels, at, height)) # combine all to df for ease of use
    
    hist_wo_legend <- lab_hist + theme_minimal() +
      geom_line(data = as.data.frame(cbind(lns, seq)),
                aes(x=seq, y = lns),
                col = "darkred",
                size = 1,
                linetype = "dashed") +
      geom_point(aes(x = predict(mod2, newdata=data.frame(z = z0)),
                     y = 0,),
                 cex = 10,
                 bg="white",
                 pch = 21) +
      geom_text(aes(x = predict(mod2, newdata=data.frame(z = z0)),
                    y = 0),
                label = "0") +
      geom_point(aes(x = a1.0,
                     y = 0),
                 cex = 10,
                 bg="white",
                 pch = 21,) +
      geom_text(aes(x = a1.0,
                    y = 0),
                label = "1") +
      geom_line(data = sun.elv.angle.df, 
                aes(x = at, 
                    y = height)) + 
      geom_label(data = sun.elv.angle.df,
                 aes(x = at, 
                     y = height, 
                     label = labels))  +
      geom_text(data = sun.elv.angle.df,
                aes(x = mean(at), 
                    y = unique(height) + 0.5 , 
                    label = "sun elevation angle (degrees)")) +
      annotate("text", 
               x = max(twl_dev)-10, #x location based on max value 
               y = unique(height) -0.2, # adjusted to be under sun elevation angle
               hjust = 0, # text is left aligned
               vjust = 1, # position regarding the other annotations
               label = paste("0. Sun elevation angle (zero):", round(90 - a1.1, 3))) + 
      annotate("text", 
               x = max(twl_dev)-10, 
               y = unique(height) - 0.25, 
               hjust = 0, 
               vjust = 3,
               label = paste("1. Sun elevation angle (median):", round(90 - a1.1, 3))) 
    
    
    if(method=="log-norm") hist_w_legend <- hist_wo_legend +
      annotate("text", x = max(twl_dev)-10, y = unique(height) - 0.3, hjust = 0, vjust = 5,
               label = paste("Log-mean:", round(fitml_ng$estimate[1], 3))) +
      annotate("text", x = max(twl_dev)-10, y = unique(height) - 0.35, hjust = 0, vjust = 7,
               label = paste("Log-sd:", round(fitml_ng$estimate[2], 3)))
    
    if(method=="gamma") hist_w_legend <- hist_wo_legend +
      annotate("text", x = max(twl_dev)-10, y = unique(height) - 0.3, hjust = 0, vjust = 5,
               label = paste("Shape:", round(fitml_ng$estimate[1], 3))) +
      annotate("text", x = max(twl_dev)-10, y = unique(height) - 0.35, hjust = 0, vjust = 7,
               label = paste("Scale:", round(fitml_ng$estimate[2], 3)))
    
    print(hist_w_legend)
  }
  
  c(a1 = as.numeric(median(z)), e0 = as.numeric(90-z0), 
    shape =  as.numeric(fitml_ng$estimate[1]), scale =  as.numeric(fitml_ng$estimate[2]))
}

##' Residency analysis using a changepoint model
##'
##' Function to discriminate between periods of residency and movement based on
##' consecutive sunrise and sunset data. The calculation is based on a
##' changepoint model (\bold{\pkg{R}} Package \code{\link{changepoint}}:
##' \code{\link{cpt.mean}}) to find multiple changepoints within the
##' data.
##'
##' The \code{cpt.mean} from the \code{R} Package \code{changepoint} is a
##' function to find a multiple changes in mean for data where no assumption is
##' made on their distribution. The value returned is the result of finding the
##' optimal location of up to Q changepoints (in this case as many as possible)
##' using the cumulative sums test statistic.
##'
##'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param quantile probability threshold for stationary site selection. Higher
##' values (above the defined quantile of all probabilities) will be considered
##' as changes in the behavior. Argmuent will only be considered if either \code{rise.prob} and/or
##' \code{set.prob} remain unspecified.
##' @param rise.prob the probability threshold for \bold{sunrise}: greater or
##' equal values indicates changes in the behaviour of the individual.
##' @param set.prob the probability threshold for \bold{sunset}: higher and
##' equal values indicates changes in the behaviour of the individual.
##' @param days a threshold for the length of stationary period. Periods smaller
##' than "days" will not be considered as a residency period
##' @param fixed ...
##' @param plot logical, if \code{TRUE} a plot will be produced
##' @param summary logical, if \code{TRUE} a summary of the results will be
##' printed
##' @return A \code{list} with probabilities for \emph{sunrise} and
##' \emph{sunset} the user settings of the probabilities and the resulting
##' stationary periods given as a \code{vector}, with the residency sites as
##' positiv numbers in ascending order (0 indicate movement/migration).
##' @note The sunrise and/or sunset times shown in the graph (if
##' \code{plot=TRUE}) represent hours of the day. However if one or both of the
##' twilight events cross midnight during the recording period the values will
##' be formed to avoid discontinuity.
##' @author Simeon Lisovski & Tamara Emmenegger
##' @seealso \code{\link{changepoint}}, \code{\link{cpt.mean}}
##' @references Taylor, Wayne A. (2000) Change-Point Analysis: A Powerful New
##' Tool For Detecting Changes.
##'
##' M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis.
##' \emph{Wiley}.
##'
##' Chen, J. and Gupta, A. K. (2000) Parametric statistical change point
##' analysis. \emph{Birkhauser}.
##' @examples
##'
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' residency <- changeLight(hoopoe2, quantile=0.9)
##'
##' @importFrom changepoint cpt.mean cpts.full pen.value.full
##' @importFrom stats na.omit quantile aggregate
##' @importFrom ggplot2 ggplot aes geom_line scale_y_continuous geom_rect theme_bw labs theme element_rect element_blank scale_x_datetime sec_axis geom_bar geom_hline
##' @importFrom patchwork wrap_plots
##' @export changeLight
changeLight <- function (tFirst, tSecond, type, twl, quantile = 0.9, rise.prob = NA, set.prob = NA, days = 5, fixed = NULL, plot = TRUE, summary = TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), 
                                                  FUN = function(x) any(class(x) != "name"))])
  
  if(is.null(fixed)) fixed <- matrix(FALSE, ncol = 2, nrow = nrow(tab))
  
  tw <- data.frame(datetime = as.POSIXct(c(as.numeric(tab$tFirst), as.numeric(tab$tSecond)), origin = "1970-01-01", "GMT"), 
                   type = c(tab$type, ifelse(tab$type == 1, 2, 1)), row = rep(1:nrow(tab), 2),
                   fixed = c(fixed[,1], fixed[,2]))
  tw <- tw[!duplicated(tw$datetime), ]
  tw <- tw[order(tw[, 1]), ]
  hours <- as.numeric(format(tw[, 1], "%H")) + as.numeric(format(tw[,1], "%M"))/60
  
  for (t in 1:2) {
    cor <- rep(NA, 24)
    for (i in 0:23) {
      cor[i + 1] <- max(abs((c(hours[tw$type == t][1], 
                               hours[tw$type == t]) + i)%%24 - (c(hours[tw$type == t], 
                                                                  hours[tw$type == t][length(hours)]) + i)%%24), na.rm = T)
    }
    hours[tw$type == t] <- (hours[tw$type == t] + (which.min(round(cor, 2))) - 1)%%24
  }
  
  sr <- tw[tw[, 2]==1, 1]
  ss <- tw[tw[, 2]==2, 1]
  rise <- hours[tw[, 2] == 1]
  set  <- hours[tw[, 2] == 2]
  CPs1 <- suppressWarnings(cpt.mean(rise, method = "BinSeg", 
                                    Q = length(rise)/2, penalty = "Manual", pen.value = 0.001, 
                                    test.stat = "CUSUM", param.estimates = FALSE))
  CPs2 <- suppressWarnings(cpt.mean(set, method = "BinSeg", 
                                    Q = length(set)/2, penalty = "Manual", pen.value = 0.001, 
                                    test.stat = "CUSUM", param.estimates = FALSE))
  N1 <- seq(1, length(rise))
  N2 <- seq(1, length(set))
  tab1 <- merge(data.frame(N = N1, prob = NA), data.frame(N = cpts.full(CPs1)[nrow(cpts.full(CPs1)),], 
                                                          prob = pen.value.full(CPs1)/2), by.x = "N", by.y = "N", all.x = T)[, -2]
  tab1[is.na(tab1[, 2]), 2] <- 0
  tab1[tw$fixed[tw$type==1],2] <- NA
  
  tab2 <- merge(data.frame(N = N2, prob = NA), data.frame(N = cpts.full(CPs2)[nrow(cpts.full(CPs2)),], 
                                                          prob = pen.value.full(CPs2)/2), by.x = "N", by.y = "N",all.x = T)[, -2]
  tab2[is.na(tab2[, 2]), 2] <- 0
  tab2[tw$fixed[tw$type==2],2] <- NA
  
  
  if (is.na(rise.prob) & is.na(set.prob)) {
    rise.prob <- as.numeric(round(as.numeric(quantile(tab1[tab1[,2] != 0, 2], probs = quantile, na.rm = TRUE)), digits = 5))
    set.prob <- as.numeric(round(as.numeric(quantile(tab2[tab2[,2] != 0, 2], probs = quantile, na.rm = TRUE)), digits = 5))
  }
  
  riseProb <- data.frame(time = tw[tw[, 2] == 1, 1], prob = tab1[,2])
  setProb  <- data.frame(time = tw[tw[, 2] == 2, 1], prob = tab2[,2])
  
  tmp02   <- data.frame(tab, rise.prob = apply(cbind(as.numeric(tab[,1]), as.numeric(tab[,2]), tab$type), 1, 
                                               function(x) ifelse(x[3]==1, riseProb$prob[which.min(abs(x[1]-as.numeric(riseProb$time)))],
                                                                  riseProb$prob[which.min(abs(x[2]-as.numeric(riseProb$time)))])),
                        set.prob = apply(cbind(as.numeric(tab[,1]), as.numeric(tab[,2]), tab$type), 1, 
                                         function(x) ifelse(x[3]==2, setProb$prob[which.min(abs(x[1]-as.numeric(setProb$time)))],
                                                            setProb$prob[which.min(abs(x[2]-as.numeric(setProb$time)))])))
  
  
  tmp02$cut   <- ifelse(apply(tmp02[, c("rise.prob", "set.prob", "type")], 1, function(x) any(ifelse(x[3]==1, x[1]>=rise.prob, x[2]>=set.prob), ifelse(x[3]==1, x[2]>=set.prob, x[1]>=rise.prob))), NA, TRUE)
  tmp02$fixed <- apply(fixed, 1, function(x) any(x))
  
  tmp02 <- cbind(tmp02, NA)
  
  s <- 1
  for (i in 1:nrow(tmp02)) {
    if(i<nrow(tmp02)) if(tmp02$fixed[i+1] & !tmp02$fixed[i]) tmp02$cut[i] <- NA
    if(i>1) if(tmp02$fixed[i-1] & !tmp02$fixed[i]) tmp02$cut[i] <- NA
    
    if(tmp02[i, 'fixed']) {
      if(i>1) if(is.na(tmp02$cut[i-1]) & !tmp02$fixed[i-1]) s <- s+1
      tmp02[i, 8] <- s} else { 
        if(i%in%c(2:(nrow(tmp02)-1))){
          if(is.na(tmp02[i - 1, 'cut']) & !is.na(tmp02[i, 'cut'])) {
            s <- s + 1
            tmp02[i, 8] <- s
          }
          if (!is.na(tmp02[i - 1, 'cut']) & !is.na(tmp02[i, 'cut'])) 
            tmp02[i, 8] <- s
        }
      }
  }
  ind01 <- aggregate(as.numeric(tmp02[!is.na(tmp02[,8]) & !tmp02$fixed,1]), by = list(tmp02[!is.na(tmp02[,8]) & !tmp02$fixed, 8]), 
                     FUN =  function(x) (x[length(x)] - x[1])/60/60/24 > days)
  tmp02[, 8] <- ifelse(tmp02[, 8]%in%c(ind01[ind01[,2],1], unique(tmp02[tmp02$fixed, 8])), tmp02[, 8], NA)
  
  s <- 1
  for (i in unique(tmp02[!is.na(tmp02[,8]),8])) {
    tmp02[!is.na(tmp02[, 8]) & tmp02[, 8] == i, 8] <- s
    s <- s + 1
  }
  
  t02 <- schedule(tmp02$tFirst, tmp02$tSecond, tmp02[,8])
  
  arr <- tmp02[!is.na(tmp02[, 8]) & !duplicated(tmp02[, 8]),]
  dep <- tmp02[!is.na(tmp02[, 8]) & !duplicated(tmp02[, 8], fromLast = T), ]
  t02$P.start <- ifelse(arr$type==1, arr$rise.prob, arr$set.prob)
  t02$P.end <- ifelse(dep$type==1, dep$rise.prob, dep$set.prob)
  t02$Days <- apply(t02, 1, function(x) round(as.numeric(difftime(x[3], x[2], units = "days")), 1))
  ds <- t02
  out <- list(riseProb = tab1[, 2], setProb = tab2[, 2], rise.prob = rise.prob, 
              set.prob = set.prob, site = ifelse(is.na(tmp02[,8]), 0 , tmp02[,8]), migTable = ds)
  
  if (plot) {
    # needed for poc_plots and sunset/sunrise, also for all xlab labels and breaks
    tw_rise <- data.frame(sr, rise)
    tw_set <- data.frame(ss, set)
    
    # Old part of code necessary for histogram + rectangles 
    mig <- out$site
    mig[mig > 0] <- 1
    
    min.r <- aggregate(as.numeric(tab$tFirst[out$site>0]),
                       by = list(out$site[out$site>0]), 
                       FUN = function(x) min(x))
    
    max.r <- aggregate(as.numeric(tab$tFirst[out$site>0]),
                       by = list(out$site[out$site>0]),
                       FUN = function(x) max(x))
    
    # New part of code necessary for histogram + rectangles
    min.r[, 2] <- as.POSIXct(min.r[, 2], 
                             origin = "1970-01-01",
                             tz = "UTC")
    
    max.r[, 2] <- as.POSIXct(max.r[, 2], 
                             origin = "1970-01-01",
                             tz = "UTC")
    
    # Histogroam without rectangles
    histogram <-  ggplot() +
      geom_line(aes(x = tab[, 1] + (tab[, 2] - tab[, 1])/2, 
                    y = ifelse(out$site > 0, 1, 0))) +
      scale_y_continuous(limits = c(0, 1.5)) 
    
    # Add rectangles to histogram
    hist_rect <- histogram + 
      geom_rect(aes(xmin = min.r[, 2],
                    xmax = max.r[, 2], 
                    ymin = 1.1, 
                    ymax = 1.4),
                fill = "grey30",
                color = "transparent",
                size = 0) +
      geom_rect(
        aes(xmin = min.r[, 2],
            xmax = max.r[, 2],
            ymin = 1.1, 
            ymax = 1.4),
        fill = ifelse(unique(out$site[out$site > 0]) %in% unique(tmp02[tmp02$fixed, 8]), "red", "transparent"),
        color = "transparent",
        size = 0,
        alpha = 0.5
      ) +
      theme_bw() +
      labs(x = "", 
           y = "") +
      theme(panel.background = element_rect(fill = "white"),
            panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_blank()) +
      scale_x_datetime(breaks = seq(min(sr), 
                                    max(sr), 
                                    length.out = 10), 
                       date_labels = "%b %y")
    
    # sunrise sunset plot 
    
    axis_increase <- 7 #increase two have sunrise and sunset on equal height in plot
    
    rise_set <- ggplot() +
      geom_line(data = tw_rise, 
                aes(x = sr, 
                    y = rise), 
                size = 2,
                color = "firebrick", 
                lwd = 0.5) +
      geom_line(data = tw_set, 
                aes(x = ss,
                    y = set + axis_increase), 
                color = "cornflowerblue", 
                size = 2, 
                lwd = 0.5) +
      scale_y_continuous(name = "Sunrise (red)",
                         sec.axis = sec_axis(~.-axis_increase, 
                                             name="Sunset (blue)")) +
      scale_x_datetime(name = "",
                       breaks = seq(min(sr), 
                                    max(sr), 
                                    length.out = 10), 
                       date_labels = "%b %y") +
      theme_bw() +
      theme(panel.background = element_rect(fill = "white"),
            panel.grid = element_blank(),
            axis.text.x = element_blank())
    
    poc_red <- ggplot() +
      geom_bar(aes(x = sr,
                   y = tab1[, 2]), 
               stat = "identity", 
               color = "firebrick") +
      labs(x = "",
           subtitle = "Sunrise", 
           y = "Probabilty of change") +
      theme_bw() +
      theme (panel.background = element_rect(fill = "white"),
             axis.text.x = element_blank(),
             panel.grid = element_blank()) +
      scale_x_datetime(breaks = seq(min(sr), 
                                    max(sr), 
                                    length.out = 10), 
                       date_labels = "%b %y") 
    
    if (is.numeric(rise.prob)) {
      poc_red +
        geom_hline(yintercept = rise.prob, 
                   lty = 2, 
                   lwd = 1.5)
    }
    
    poc_blue <- ggplot()+
      geom_bar(aes(x = ss,
                   y = tab2[, 2]),
               stat = "identity",
               color = "cornflowerblue") +
      labs(x = "Time",
           subtitle = "Sunset",
           y = "Probabilty of change") + 
      theme_bw() +
      theme (panel.background = element_rect(fill = "white"),
             panel.grid = element_blank()) +
      scale_x_datetime(breaks = seq(min(sr), 
                                    max(sr), 
                                    length.out = 10), 
                       date_labels = "%b %y") 
    
    if (is.numeric(set.prob)) 
      poc_blue +
      geom_hline(yintercept = set.prob, lty = 2, lwd = 1.5)
    
    print(wrap_plots(hist_rect / (rise_set) / (poc_red) / (poc_blue))) # use of patchwork syntax to show all plots below each other
  }
  
  if (summary) {
    i.sum.Cl(out)
  }
  return(out)
}

#' siteEstimate
#'
#' ...
#'
#' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
#' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
#' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
#' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type}
#' @param degElevation the sun elevation angle (in degrees) that defines twilight (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector}.
#' @param method \code{character} string; only \code{gamma} and \code{log-normal} are implemented.
#' @param parms a \code{vector} describing the two parameters of the error density distribution (defined by \code{method}).
#' @param xlim the longitudinal boundaries for which the likelihood will be calculated.
#' @param ylim the latitudinal boundaries for which the likelihood will be calculated.
#' @param res the spatial resolution in degrees.
#' @return A \code{list} with ...
#' @author Simeon Lisovski
#'
#' @export siteEstimate
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ parRapply stopCluster detectCores
siteEstimate <- function(tFirst, tSecond, type, twl, 
                         degElevation, 
                         method = "gamma", parms = c(3.3, 0.8), 
                         xlim = c(-180, 180), 
                         ylim = c(-90, 90), res = c(0.5, 0.5)) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  tw <- data.frame(Twilight = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   Rise = c(ifelse(tab$type==1, TRUE, FALSE), ifelse(tab$type == 1, FALSE, TRUE)))
  tw <- tw[!duplicated(tw$Twilight),]
  tw <- tw[order(tw[,1]),]
  
  
  loglik <- function(crds, Twilight, Rise, degElevation, method, parms) {
    t.tw <- twilight(Twilight, lon = crds[1], lat = crds[2], 
                     rise = ifelse(Rise, TRUE, FALSE), zenith = 90-degElevation, 
                     iters = 6)
    diff.sr <- as.numeric(difftime(Twilight[Rise], t.tw[Rise], units = "mins"))
    diff.ss <- as.numeric(difftime(t.tw[!Rise], Twilight[!Rise], units = "mins"))
    if(method=="log-norm") {
      return(-sum(dlnorm(c(diff.sr, diff.ss), parms[1], parms[2], log = T), na.rm = T))
    }
    if(method=="gamma") {
      return(-sum(dgamma(c(diff.sr, diff.ss), parms[1], parms[2], log = T), na.rm = T))
    }
  }
  
  
  lon <- seq(xlim[1], xlim[2], by = res[1])
  lat <- rev(seq(ylim[1], ylim[2], by = res[2]))
  
  crdsm <- data.frame(lon = rep(lon, length(lat)), 
                      lat = rep(lat, each = length(lon)))
  
  out_A  <- array(dim = c(length(lat), length(lon), length(degElevation)))
  colnames(out_A) <- lon
  rownames(out_A) <- lat
  out_ML <- matrix(NA, ncol = 2, nrow = length(degElevation))
  
  
  eq.ind <- which(format(tw$Twilight, "%m/%d")%in%c("03/21", "09/21"))
  
  if(length(eq.ind)>0) {
    twl.sp <- split(tw, f = ifelse(1:nrow(tw)<min(eq.ind), 1, ifelse(1:nrow(tw)>max(eq.ind), 2, NA)))
  } else {
    twl.sp <- list(tw)
  }
  
  mycl <- makeCluster(detectCores()-1)
  tmp  <- clusterSetRNGStream(mycl)
  tmp  <- clusterExport(mycl, c("loglik"), envir=environment())
  tmp  <- clusterEvalQ(mycl, library("GeoLight")) 
  
  for(i in 1:length(degElevation)) {
    nll0  <- lapply(twl.sp, function(x) parRapply(mycl, crdsm, FUN = loglik, Twilight = x$Twilight, Rise = x$Rise, 
                                                  degElevation = degElevation[i], method = method, parms = parms))
    nll   <- apply(do.call("cbind", nll0), 1, function(x) ifelse(all(is.infinite(abs(x))), Inf, sum(x)))
    
    out_A[,,i] <- matrix(suppressWarnings((max(nll[is.finite(nll)])-nll)/sum(nll[is.finite(nll)], na.rm = T)), 
                         ncol = length(lon), nrow = length(lat), byrow = T)
    if(any(is.finite(abs(out_A[,,i])))) {
      out_ML[i,] <- cbind(as.numeric(colnames(out_A)[which(out_A[,,i] == max(out_A[,,i], na.rm = T), arr.ind = TRUE)[2]]),
                          as.numeric(rownames(out_A)[which(out_A[,,i] == max(out_A[,,i], na.rm = T), arr.ind = TRUE)[1]]))
    }
  }
  
  stopCluster(mycl)
  
  list(SunElevation = degElevation,
       Estimate = out_A,
       mlLoc = out_ML)
}

#' Function to merge sites
#'
#' The \code{\link{changeLight}} functions provides a vector grouping the twilight times
#' into stationary (>0) and movement (0) periods. This function was written to enable the user
#' to merge sites based on the distance between consequtive sites. NOTE: The function requires
#' position estimate and desicison on whether sites should be merged will be made based on
#' the defined \code{distance}, the \code{cutoff} values and the provided positions. The analysis
#' is this dependent on the accuracy of the position estiamtes and should be applied to positons that
#' were estimated using a sensible sun elevation angle.
#'
#' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
#' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
#' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
#' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type}
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0. This \code{vector} will be used as the initial state.
#' @param degElevation the sun elevation angle (in degrees) that defines twilight (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst} or \code{nrow(x)}.
#' @param distThreshold a \code{numerical} value defining the threshold of the distance under 
#' which consequtive sites should be merged (in km).
#' @param fixed ...
#' @param alpha mean and standard variation for position optimization process.
#' @param plot \code{logical}, if TRUE a plot comparing the inital and the final site selection.
#' @return A \code{vector} with the merged site numbers
#' @author Simeon Lisovski
#'
#' @export mergeSites
#' @importFrom fields rdist.earth
#' @importFrom graphics abline axis lines mtext par plot points rect 
#' @importFrom stats optim dnorm
mergeSites <- function(tFirst, tSecond, type, twl, site, degElevation, distThreshold = 250, 
                       fixed = NULL, alpha = c(0, 15), plot = TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), 
                                                  FUN = function(x) any(class(x) != "name"))])
  
  if(is.null(fixed)) fixed <- matrix(FALSE, ncol = 2, nrow = nrow(tab))
  fixed.ind <- apply(fixed, 1, function(x) any(x))
  
  site0 <- site
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type = c(tab$type, ifelse(tab$type == 1, 2, 1)))
  tw <- tw[!duplicated(tw$datetime), ]
  tw <- tw[order(tw[, 1]), ]
  tw <- tw[1:nrow(tab),]
  
  crds0 <- coord(tab, degElevation = degElevation, note = F)
  tab$lon <- crds0[,1]
  tab$lat <- crds0[,2]
  tw$lon <- crds0[,1]
  tw$lat <- crds0[,2]
  
  lonlim <- range(crds0[, 1], na.rm = T)
  lon.seq <- seq(lonlim[1] - 1, lonlim[2] + 1, by = 1)
  latlim <- range(crds0[, 2], na.rm = T)
  lat.seq <- seq(latlim[1] - 1, latlim[2] + 1, by = 1)
  mod <- function(x) {
    loglik <- function(crds) {
      t.tw <- twilight(x$tFirst, lon = crds[1], lat = crds[2], 
                       rise = ifelse(x$type == 1, TRUE, FALSE), zenith = 90 - 
                         degElevation, iters = 6)
      diff <- as.numeric(difftime(x$tFirst, t.tw, units = "mins"))
      -sum(dnorm(diff, alpha[1], alpha[2], log = T), na.rm = T)
    }
    fit0 <- optim(cbind(median(x$lon, na.rm = T), median(x$lat, na.rm = T)), loglik, lower = cbind(lonlim[1], 
                                                                                                   latlim[1]), upper = cbind(lonlim[2], latlim[2]), 
                  method = "L-BFGS-B", hessian = T)
    fit  <- optim(cbind(fit0$par[1], fit0$par[2]), loglik, lower = cbind(lonlim[1], 
                                                                         latlim[1]), upper = cbind(lonlim[2], latlim[2]), 
                  method = "L-BFGS-B", hessian = T)
    fisher_info <- solve(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    lon.lower <- c(fit$par[1] - 1.96 * prop_sigma)[1]
    lat.lower <- c(fit$par[2] - 1.96 * prop_sigma)[4]
    lon.upper <- c(fit$par[1] + 1.96 * prop_sigma)[1]
    lat.upper <- c(fit$par[2] + 1.96 * prop_sigma)[4]
    matrix(c(fit$par[1], fit$par[2], lon.lower, lat.lower, 
             lon.upper, lat.upper), ncol = 6)
  }
  out <- data.frame(site = unique(site[site != 0]), t(sapply(split(tab[site != 0, ], f = site[site != 0]), mod)))
  rep = TRUE
  ite = 1
  repeat {
    for (i in site[site != 0 & !duplicated(site) & !fixed.ind & !site%in%(unique(site[fixed.ind])-1)]) {
      if(i==max(out$site)) break
      dist0 <- rdist.earth(out[which(out[,1]==i):(which(out[,1]==i) + 1), 2:3])[2, 1]
      if (dist0 <= distThreshold) 
        break
    }
    if (i < max(site[site != 0 & !duplicated(site) & !fixed.ind & !site%in%(unique(site[fixed.ind])-1)])) {
      site[(which(site == i)[1]):(which(site == (i + 1))[sum(site == (i + 1))])] <- i
      site[which(site > i)] <- site[which(site > i)] - 1
    } else rep = FALSE
    out <- data.frame(site = unique(site[site != 0 & !fixed.ind]), t(sapply(split(tab[site != 0 & !fixed.ind, ], 
                                                                                  f = site[site != 0 & !fixed.ind]), mod)))
    if (!rep) 
      break
    else ite <- ite + 1
  }
  
  if(any(fixed)) {
    fs <- site[site != 0 & !duplicated(site) & fixed.ind & !site%in%(unique(site[fixed.ind])-1)]
    out.temp <- as.data.frame(cbind(fs, matrix(NA, nrow = length(fs), ncol = ncol(out)-1)))
    names(out.temp) <- names(out)
    out <- rbind(out, out.temp)
    out <- out[order(out[,1]),]
  }
  
  
  if (plot) {
    hours0 <- as.numeric(format(tw[, 1], "%H")) + as.numeric(format(tw[, 1], "%M"))/60
    crd0 <- out[match(site, out$site), 2:3]
    crd0[!is.na(crd0[,1]),] <- crds0[!is.na(crd0[,1]),]
    hours1 <- twilight(tw[, 1], rise = ifelse(tw[, 2] == 1, TRUE, FALSE), zenith = 90 - degElevation, 
                       lon = out[match(site, out$site), 2], lat = out[match(site, out$site), 3])
    hours1 <- as.numeric(format(hours1, "%H")) + as.numeric(format(hours1, "%M"))/60
    hours2 <- twilight(tw[, 1], rise = ifelse(tw[, 2] == 1, TRUE, FALSE), zenith = 90 - degElevation, 
                       lon = out[match(site, out$site), 4], lat = out[match(site, out$site), 3])
    
    hours2 <- as.numeric(format(hours2, "%H")) + as.numeric(format(hours2,"%M"))/60
    hours3 <- twilight(tw[, 1], rise = ifelse(tw[, 2] == 1, TRUE, FALSE), zenith = 90 - degElevation, 
                       lon = out[match(site, out$site), 6], lat = out[match(site, out$site), 3])
    hours3 <- as.numeric(format(hours3, "%H")) + as.numeric(format(hours3,"%M"))/60
    
    
    for (t in 1:2) {
      cor <- rep(NA, 24)
      for (i in 0:23) {
        cor[i + 1] <- max(abs((c(hours0[tw$type == t][1], 
                                 hours0[tw$type == t]) + i)%%24 - 
                                (c(hours0[tw$type == t], 
                                   hours0[tw$type == t][length(hours0)]) +i)%%24), 
                          na.rm = T)
      }
      hours0[tw$type == t] <- (hours0[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
      hours1[tw$type == t] <- (hours1[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
      hours2[tw$type == t] <- (hours2[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
      hours3[tw$type == t] <- (hours3[tw$type == t] + (which.min(round(cor,2))) - 1)%%24
    }
    
    opar <- par(mfrow = c(5, 1), oma = c(5, 0, 0, 0), mar = c(1.5,5, 1, 1))
    mig1 <- site0
    mig1[mig1 > 0] <- 1
    mig2 <- site
    mig2[mig2 > 0] <- 1
    plot(tw[, 1], ifelse(mig2 > 0, 1, 0), type = "l", yaxt = "n", 
         ylab = NA, ylim = c(0, 1.5), col = "firebrick", lwd = 2, 
         xaxt = "n")
    lines(tw[, 1], ifelse(mig1 > 0, 1, 0), type = "l", lty = 2)
    rect(tw[site > 0 & !duplicated(site), 1], 1.1, 
         tw[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, col = "grey90", lwd = 0)
    rect(tw[site > 0 & !duplicated(site), 1], 1.1, 
         tw[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, 
         col = ifelse(apply(fixed[site>0,], 1, function(x) any(x))[!duplicated(site[site>0])], "red", "transparent"), 
         density = 60)
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[tw[, 2] == 1, 1], hours1[tw[, 2] == 1], type = "l", 
         lwd = 2, col = "firebrick", ylab = "Sunrise (red)", 
         xlim = range(tw[, 1]), ylim = range(hours0[tw[, 2] == 1]), xaxt = "n")
    lines(tw[tw[, 2] == 1, 1], hours2[tw[, 2] == 1], type = "l", 
          lwd = 1, lty = 2)
    lines(tw[tw[, 2] == 1, 1], hours3[tw[, 2] == 1], type = "l", 
          lwd = 1, lty = 2)
    points(tw[tw[, 2] == 1 & !fixed.ind, 1], hours0[tw[, 2] == 1 & !fixed.ind], cex = 0.5, 
           pch = 21, col = "black", bg = "firebrick", lwd = 0.5)
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[tw[, 2] == 2, 1], hours1[tw[, 2] == 2], type = "l", 
         lwd = 2, col = "cornflowerblue", ylab = "Sunset (blue)", 
         xlim = range(tw[, 1]), ylim = range(hours0[tw[, 2] == 
                                                      2]), xaxt = "n")
    lines(tw[tw[, 2] == 2, 1], hours2[tw[, 2] == 2], type = "l", 
          lwd = 1, lty = 2)
    lines(tw[tw[, 2] == 2, 1], hours3[tw[, 2] == 2], type = "l", 
          lwd = 1, lty = 2)
    points(tw[tw[, 2] == 2 & !fixed.ind, 1], hours0[tw[, 2] == 2 & !fixed.ind], cex = 0.5, 
           pch = 21, col = "black", bg = "cornflowerblue", lwd = 0.5)
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[, 1], ifelse(fixed.ind, NA, crds0[, 1]), type = "o", pch = 16, cex = 0.5, 
         xaxt = "n", ylab = "Longitude", cex.lab = 1.7, xlab = "")
    abline(v = c(tw[site0 > 0 & !duplicated(site0), 1], 
                 tw[site0 > 0 & !duplicated(site0, fromLast = T), 1]), lty = 2)
    abline(v = c(tw[site > 0 & !duplicated(site), 1], 
                 tw[site > 0 & !duplicated(site, fromLast = T), 1]), lwd = 1.5, 
           col = "firebrick")
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = FALSE)
    plot(tw[, 1], ifelse(fixed.ind, NA, crds0[, 2]), type = "o", pch = 16, cex = 0.5, 
         xaxt = "n", ylab = "Latitude", cex.lab = 1.7, xlab = "")
    abline(v = c(tw[site0 > 0 & !duplicated(site0), 1], 
                 tw[site0 > 0 & !duplicated(site0, fromLast = T), 1]), lty = 2)
    abline(v = c(tw[site > 0 & !duplicated(site), 1], 
                 tw[site > 0 & !duplicated(site, fromLast = T), 1]), lwd = 1.5, 
           col = "firebrick")
    axis(1, at = seq(tw[1, 1], tw[nrow(tw), 1], length = 10), 
         labels = format(seq(tw[1, 1], tw[nrow(tw), 1], length = 10), "%d-%b"))
    mtext("Date", 1, outer = T, line = 1.6, cex = 1.2)
    par(opar)
  }
  names(out) <- c("site", "Lon", "Lat", "Lon.upper", "Lat.upper", 
                  "Lon.lower", "Lat.lower")
  list(site = site, summary = out)
}



#' Imporved mergeSites function
#'
#' The \code{\link{changeLight}} functions provides a vector grouping the twilight times
#' into stationary (>0) and movement (0) periods. The mergeSites functions allows to merge
#' these stationary periods in case some consecutive sites were separated by e.g. outliers
#' or strong shading events. In a nutshell, the function estimates a likelihood surface for each
#' stationary period that is based on the siteEstimate principle (e.g. fitting sunrise and sunset times
#' of locations to the recorded sunrise and sunsets and finding the coordinates that results in
#' the best fit).
#' The decision on whether two consecutive sites should be merged together is made by the provided
#' `distThreshold` and the overlap of the 95% credible intervals of the location (e.g. the distance
#' between two consecutive sites needs to be smaller than the `distThreshold` AND the median locations need
#' to overlap with the credible intervals of the two location estimates to trigger merging of the sites).
#' The function also requires a sun elevation angle. HOWEVER, the angle defining the zero in the
#' log-normal distribution of the twilight error is needed and not the sun elevation for the median
#' twilight error which is needed for the location estimation using the function `coord`. The function `getElevation`
#' provides both, the median and the zero sun elevation angle as well as the parameters for the log-normal
#' distribution that also required in the `mergeSites2` function.
#'
#' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
#' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
#' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
#' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type}
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0. This \code{vector} will be used as the initial state.
#' @param degElevation the sun elevation angle (in degrees) that defines twilight (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst} or \code{nrow(x)}.
#' @param distThreshold a \code{numerical} value defining the threshold of the distance under 
#' which consecutive sites should be merged (in km).
#' 
#' @param fixed a logical vector indicating locations that should not be merged.
#' @param alpha log-mean and log-sd of the twilight error (see: \code{\link{getElevation}}).
#' @param method either, "gamma" or "log-normal". Depending on the error distribution used in \code{\link{getElevation}}.
#' @param map A 'SpatialPolygonDataFrame' of the world. E.g., use 'wrld_simpl' from the maptools package.
#' @param mask Either 'land' or 'ocean'.
#' @param plot \code{logical}, if TRUE a plot comparing the initial and the final site selection.
#' @param ggplot logical. if TRUE plots will be created using ggplot2
#' @return A \code{vector} with the merged site numbers
#' @author Simeon Lisovski
#'
#' @export mergeSites2
#' @importFrom fields rdist.earth
#' @importFrom terra rast rasterize crds xFromCell yFromCell ncell
#' @importFrom graphics abline axis lines mtext par plot points rect 
#' @importFrom stats optim dnorm
#' @importFrom utils data 
#' @importFrom parallel makeCluster clusterSetRNGStream parRapply stopCluster detectCores
#' @importFrom ggplot2 scale_x_datetime theme_bw theme element_blank labs geom_point geom_vline
#' @importFrom patchwork wrap_plots
mergeSites2 <- function(tFirst, tSecond, type, twl, site, degElevation,
                        distThreshold = 250, fixed = NULL, alpha = c(2.5, 1), method = "gamma", map, mask = "land", plot = TRUE, ggplot = TRUE) {
  
  if(!(method%in%c("gamma", "log-norm"))) stop("Method can only be `gamma` or `log-norm`.")
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), 
                                                  FUN = function(x) any(class(x) != "name"))])
  
  mycl <- makeCluster(detectCores()-1)
  tmp  <- clusterSetRNGStream(mycl)
  
  if(is.null(fixed)) fixed <- matrix(FALSE, ncol = 2, nrow = nrow(tab))
  if(is.null(fixed)) stop("fixed is NULL")
  
  fixed.ind <- apply(fixed, 1, function(x) any(x))
  
  site0 <- site
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type     = c(tab$type, ifelse(tab$type == 1, 2, 1)), site = site, fixed = fixed.ind)
  tw <- tw[!duplicated(tw$datetime), ]
  tw <- tw[order(tw[, 1]), ]
  tw <- tw[1:nrow(tab),]
  tw$Rise <- ifelse(tw$type==1, TRUE, FALSE)
  site.old <- tw$site  
  
  crds0 <- coord(tab, degElevation = degElevation, note = F)
  tab$lon <- crds0[,1]
  tab$lat <- crds0[,2]
  tw$lon  <- crds0[,1]
  tw$lat  <- crds0[,2]
  
  lonlim <- range(crds0[, 1], na.rm = T)
  lon.seq <- seq(lonlim[1] - 1, lonlim[2] + 1, by = 1)
  latlim <- range(crds0[, 2], na.rm = T)
  lat.seq <- seq(latlim[1] - 1, latlim[2] + 1, by = 1)
  
  tmp  <- clusterEvalQ(mycl, {
    library(GeoLight) 
  })
  
  mod <- function(x) {
    
    lons <- seq(-180, 180, by = 2.5)
    lats <- seq(-75, 75, by = 2.5)
    
    crdsT <- expand.grid(lons, lats)
    
    ll    <- parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE), degElevation = degElevation, parms = alpha, method = method)
    ll    <- ll/max(ll, na.rm = T)
    
    r0     <- rast(xmin =-180, xmax = 180, ymin = -75, ymax = 75, nrows = length(lats), ncols = length(lons))
    r      <- rasterize(as.matrix(crdsT[!is.na(ll),]), r0, ll[!is.na(ll)])
    # plot(r, breaks = seq(0.4, 1, length = 100), col = rev(terrain.colors(99)))
    # points(xTab[[1]]$lon, xTab[[1]]$lat)
    # points(lon.calib, lat.calib, pch = 21, cex = 2, bg = "white")
    
    crdsT <- crds(r)[r[]>0.4,]           
    r0     <- rast(xmin = min(crdsT[,1])-3,
                   xmax = max(crdsT[,1])+3, 
                   ymin = min(crdsT[,2])-3, 
                   ymax = max(crdsT[,2])+3, 
                   res = 0.5)
    if(!is.null(mask)) {
      maskR  <- rasterize(map, r0)
      if(mask=="land")  crdsT  <- crds(r0)[!is.na(maskR[]),]
      
      if(mask=="ocean") crdsT  <- crds(r0)[is.na(maskR[]),]
    } else {
      crdsT  <- crds(r0)
    }
    
    ll.sr  <- parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE), degElevation = degElevation, parms = alpha, method = method, twilight = "sr")
    ll.ss    <- parRapply(mycl, crdsT, FUN = gl.loglik, Twilight = x$datetime, Rise = ifelse(x$type==1, TRUE, FALSE), degElevation = degElevation, parms = alpha, method = method, twilight = "ss")
    ll <- apply(cbind(ll.sr/max(ll.sr, na.rm = T), ll.ss/max(ll.sr, na.rm = T)), 1, function(x) ifelse(any(x<=0.0000001), NA, sum(x)))
    centre <- crdsT[which.max(ll),]
    
    r      <- rasterize(as.matrix(crdsT[!is.na(ll),]), r0, ll[!is.na(ll)])
    r[] <- r[]/max(r[], na.rm = T)
    # plot(r)
    # contour(r, add = T, levels = c(1, 0.495))
    
    crdsRange1.x <- xFromCell(r, 1:ncell(r))[!is.na(r[]) & r[]>0.495]
    crdsRange1.y <- yFromCell(r, 1:ncell(r))[!is.na(r[]) & r[]>0.495]
    crdsRange1 <- cbind(crdsRange1.x, crdsRange1.y)
    
    crdsRange2.x <- xFromCell(r, 1:ncell(r))[!is.na(r[]) & r[]>0.7]
    crdsRange2.y <- yFromCell(r, 1:ncell(r))[!is.na(r[]) & r[]>0.7]
    crdsRange2 <- cbind(crdsRange2.x, crdsRange2.y)
    
    matrix(c(centre[1], centre[2], min(crdsRange1[,1]), min(crdsRange2[,1]), max(crdsRange2[,1]), max(crdsRange1[,1]),
             min(crdsRange1[,2]), min(crdsRange2[,2]), max(crdsRange2[,2]), max(crdsRange1[,2])), ncol = 10)
    
  }
  
  tw$merge = FALSE
  xTab <- split(tw, f = tw$site)
  x0 <- xTab[[1]]
  xTab <- xTab[-1]
  
  sm <- matrix(ncol = 10, nrow = max(site.old))
  repeat{
    
    for(i in 1:(length(xTab)-1)) {
      
      if(all(!xTab[[i]]$merge)) {
        
        cat(sprintf('\n comparing site %d and %d', i, i+1)) 
        out  <- do.call("rbind", lapply(xTab[c(i, i+1)], function(x) mod(x)))
        sm[i,] <- out[1,]
        
        # plot(out[,1], out[,2], xlim = range(out[,c(1,3,4)]), ylim = range(out[,c(2,5,6)]))
        # plot(wrld_simpl, add = T)
        # arrows(out[,1], out[,5], out[,1], out[,6], length = 0)
        # arrows(out[,3], out[,2], out[,4], out[,2], length = 0)
        
        if(all(!xTab[[i]]$fixed)) {
          
          dist <- rdist.earth(out[,1:2])[2,1]
          
          if(dist<distThreshold & 
             (out[2,3] < out[1,1] & out[2,6] > out[1,1]) & (out[2,1] > out[1,3] & out[2,1] < out[1,6]) &
             (out[2,7] < out[1,2] & out[2,10] > out[1,2]) & (out[2,2] > out[1,7] & out[2,2] < out[1,10])) {
            
            tw$site[min(which(tw$site==i)):max(which(tw$site==(i+1)))] <- i
            tw$site[tw$site>0 & tw$site>i] <-  tw$site[tw$site>0 & tw$site>i]-1
            
            xTab <- split(tw, f = tw$site)
            x0 <- xTab[[1]]
            xTab <- xTab[-1]
            
            cat(".... merge (back to site 1)")
            
            break
            
          } else {
            tw$merge[tw$site==i] <- TRUE
            cat(".... no action")
          }   
          
        }
      }
    }
    if(i==(length(xTab)-1)) {
      xTab <- split(tw, f = tw$site)
      break
    }
  }
  
  stopCluster(mycl)
  
  sm  <- cbind(1:sum(!is.na(sm[,1])), sm[!is.na(sm[,1]),])
  
  out <- do.call("rbind", xTab)[,c(1,5,3,4)]
  out <- out[order(out$datetime),]       
  
  
  if (plot) {
    site   <- out$site
    hours0 <- as.numeric(format(out[, 1], "%H")) + as.numeric(format(out[, 1], "%M"))/60
    crd0   <- sm[match(out$site, sm[,1]), 2:3]
    # crd0[is.na(crd0[,1]),] <- crds0[is.na(crd0[,1]),]
    hours1 <- twilight(out$datetime, rise = out$Rise, zenith = 90 - degElevation, 
                       lon = crd0[,1], lat = crd0[,2])
    hours1 <- as.numeric(format(hours1, "%H")) + as.numeric(format(hours1, "%M"))/60
    hours2 <- twilight(out$datetime, rise = out$Rise, zenith = 90 - degElevation, 
                       lon = crd0[,1], lat = crd0[,2])
    
    hours2 <- as.numeric(format(hours2, "%H")) + as.numeric(format(hours2,"%M"))/60
    hours3 <- twilight(out$datetime, rise = out$Rise, zenith = 90 - degElevation, 
                       lon = crd0[,1], lat = crd0[,2])
    hours3 <- as.numeric(format(hours3, "%H")) + as.numeric(format(hours3,"%M"))/60
    
    
    for (t in c(TRUE, FALSE)) {
      cor <- rep(NA, 24)
      for (i in 0:23) {
        cor[i + 1] <- max(abs((c(hours0[out$Rise==t][1], 
                                 hours0[out$Rise == t]) + i)%%24 - 
                                (c(hours0[out$Rise == t], 
                                   hours0[out$Rise == t][length(hours0)]) +i)%%24), 
                          na.rm = T)
      }
      hours0[out$Rise == t] <- (hours0[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
      hours1[out$Rise == t] <- (hours1[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
      hours2[out$Rise == t] <- (hours2[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
      hours3[out$Rise == t] <- (hours3[out$Rise == t] + (which.min(round(cor,2))) - 1)%%24
    }
    
    
    opar <- par(mfrow = c(5, 1), oma = c(5, 0, 0, 0), mar = c(1.5,5, 1, 1))
    mig1 <- site.old
    mig1[mig1 > 0] <- 1
    mig2 <- out$site
    mig2[mig2 > 0] <- 1
    
    if(ggplot) {
      hist_rec <- ggplot() +
        geom_line(aes(x = out[, 1], y = ifelse(mig2 > 0, 1, 0)), 
                  col = "firebrick", linewidth = 1) +
        geom_line(aes(x = out[, 1], y = ifelse(mig1 > 0, 1, 0)), 
                  linetype = 2) +
        geom_rect(
          aes(xmin = out[site > 0 & !duplicated(site), 1], xmax = out[site > 0 & !duplicated(site, fromLast = T), 1], ymin = 1.1, ymax = 1.4),
          fill = "grey30",
          color = "transparent",
          size = 0
        ) +
        geom_rect(
          aes(xmin = out[site > 0 & !duplicated(site), 1], xmax = out[site > 0 & !duplicated(site, fromLast = T), 1], ymin = 1.1, ymax = 1.4),
          color = ifelse(apply(fixed[site>0, ], 1, function(x) any(x))[!duplicated(site[site>0])], "red", "transparent"), 
          alpha = 0.6) +
        scale_x_datetime(breaks = seq(out[1, 1], out[nrow(out), 1], length.out = 10), date_labels = "%b %d") +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        ) 
      
      red_ts <- ggplot() +
        geom_line(aes(x = out[out$Rise, 1], y = hours1[out$Rise]), 
                  col = "firebrick", linewidth = 0.5, na.rm = TRUE) +
        geom_line(aes(x = out[out$Rise, 1], y = hours2[out$Rise]), 
                  linetype = 2, linewidth = 0.5, na.rm = TRUE) +
        geom_line(aes(x = out[out[, 2] == 1, 1], y = hours3[out[, 2] == 1]), 
                  linetype = 2, linewidth = 0.5, na.rm = TRUE) +
        labs(y = "Sunrise (red)") +
        geom_point(aes(out[out[, 2] == 1 & !fixed.ind, 1], hours0[out[, 2] == 1 & !fixed.ind]), color = "black", fill = "firebrick", shape = 21, size = 0.5) +
        scale_x_datetime(breaks = seq(out[1, 1], out[nrow(out), 1], length.out = 10), date_labels = "%b %d") +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        ) 
      
      blue_ts <- ggplot() +
        geom_line(aes(x = out[!out$Rise, 1], y = hours1[!out$Rise]), 
                  col = "cornflowerblue", linewidth = 0.5, na.rm = TRUE) +
        geom_line(aes(x = out[!out$Rise, 1], y = hours2[!out$Rise]), 
                  linetype = 2, linewidth = 0.5, na.rm = TRUE) +
        geom_line(aes(x = out[!out[, 2] == 1, 1], y = hours3[!out[, 2] == 1]), 
                  linetype = 2, linewidth = 0.5, na.rm = TRUE) +
        geom_point(aes(out[!out[, 2] == 1 & !fixed.ind, 1], hours0[!out[, 2] == 1 & !fixed.ind]), color = "black", fill = "cornflowerblue", shape = 21, size = 0.5) +
        labs(y = "Sunset (blue)", x = "") +
        scale_x_datetime(breaks = seq(out[1, 1], out[nrow(out), 1], length.out = 10), date_labels = "%b %d") + 
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        ) 
      
      
      lon_plot <- ggplot() +
        geom_point(aes(out[, 1,], ifelse(fixed.ind, NA, crds0[, 1])), size = 0.5, na.rm = TRUE) +
        geom_line(aes(out[, 1,], ifelse(fixed.ind, NA, crds0[, 1]))) + 
        geom_vline(
          xintercept = c(
            out[site.old > 0 & !duplicated(site.old), 1],
            out[site.old > 0 & !duplicated(site.old, fromLast = TRUE), 1]
          ),
          linetype = "dashed"
        ) +
        geom_vline(xintercept = c(
          out[out$site  > 0 & !duplicated(out$site), 1],
          out[out$site  > 0 & !duplicated(out$site, fromLast = T), 1]),
          col = "firebrick") +
        scale_x_datetime(breaks = seq(out[1, 1], out[nrow(out), 1], length.out = 10), date_labels = "%b %d") + 
        labs(y = "Longitude", x = "") +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank()
        ) 
      
      lat_plot <- ggplot() +
        geom_point(aes(out[, 1,], ifelse(fixed.ind, NA, crds0[, 2])), size = 0.5, na.rm = TRUE) +
        geom_line(aes(out[, 1,], ifelse(fixed.ind, NA, crds0[, 2]))) + 
        geom_vline(
          xintercept = c(
            out[site.old > 0 & !duplicated(site.old), 1],
            out[site.old > 0 & !duplicated(site.old, fromLast = TRUE), 1]
          ),
          linetype = "dashed"
        ) +
        geom_vline(xintercept = c(
          out[out$site  > 0 & !duplicated(out$site), 1],
          out[out$site  > 0 & !duplicated(out$site, fromLast = T), 1]),
          col = "firebrick") +
        scale_x_datetime(breaks = seq(out[1, 1], out[nrow(out), 1], length.out = 10), date_labels = "%b %d") + 
        labs(y = "Lattidue", x = "") +
        theme_bw()
      
      print(wrap_plots( hist_rec / red_ts / blue_ts / lon_plot / lat_plot ))
    
    } else {
      plot(out[, 1], ifelse(mig2 > 0, 1, 0), type = "l", yaxt = "n", 
           ylab = NA, ylim = c(0, 1.5), col = "firebrick", lwd = 2, 
           xaxt = "n")
      lines(out[, 1], ifelse(mig1 > 0, 1, 0), type = "l", lty = 2)
      rect(out[site > 0 & !duplicated(site), 1], 1.1, 
           out[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, col = "grey90", lwd = 0)
      rect(out[site > 0 & !duplicated(site), 1], 1.1, 
           out[site > 0 & !duplicated(site, fromLast = T), 1], 1.4, 
           col = ifelse(apply(fixed[site>0,], 1, function(x) any(x))[!duplicated(site[site>0])], "red", "transparent"), 
           density = 60)
      axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
           labels = FALSE)
      
      plot(out[out$Rise, 1], hours1[out$Rise], type = "l", 
           lwd = 2, col = "firebrick", ylab = "Sunrise (red)", 
           xlim = range(out[, 1]), ylim = range(hours0[out[, 2] == 1]), xaxt = "n")
      lines(out[out$Rise, 1], hours2[out$Rise], type = "l", 
            lwd = 1, lty = 2)
      lines(out[out[, 2] == 1, 1], hours3[out[, 2] == 1], type = "l", 
            lwd = 1, lty = 2)
      points(out[out[, 2] == 1 & !fixed.ind, 1], hours0[out[, 2] == 1 & !fixed.ind], cex = 0.5, 
             pch = 21, col = "black", bg = "firebrick", lwd = 0.5)
      axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
           labels = FALSE)
      plot(out[!out$Rise, 1], hours1[!out$Rise], type = "l", 
           lwd = 2, col = "cornflowerblue", ylab = "Sunset (blue)", 
           xlim = range(out[, 1]), ylim = range(hours0[!out$Rise]), xaxt = "n")
      lines(out[!out$Rise, 1], hours2[!out$Rise], type = "l", 
            lwd = 1, lty = 2)
      lines(out[!out$Rise, 1], hours3[!out$Rise], type = "l", 
            lwd = 1, lty = 2)
      points(out[!out$Rise & !fixed.ind, 1], hours0[!out$Rise & !fixed.ind], cex = 0.5, 
             pch = 21, col = "black", bg = "cornflowerblue", lwd = 0.5)
      axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
           labels = FALSE)
      
      plot(out[, 1], ifelse(fixed.ind, NA, crds0[, 1]), type = "o", pch = 16, cex = 0.5, 
           xaxt = "n", ylab = "Longitude", cex.lab = 1.7, xlab = "")
      abline(v = c(out[site.old > 0 & !duplicated(site.old), 1], 
                   out[site.old > 0 & !duplicated(site.old, fromLast = T), 1]), lty = 2)
      abline(v = c(out[out$site  > 0 & !duplicated(out$site), 1], 
                   out[out$site  > 0 & !duplicated(out$site, fromLast = T), 1]), lwd = 1.5, col = "firebrick")
      axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
           labels = FALSE)
      plot(out[, 1], ifelse(fixed.ind, NA, crds0[, 2]), type = "o", pch = 16, cex = 0.5, 
           xaxt = "n", ylab = "Latitude", cex.lab = 1.7, xlab = "")
      abline(v = c(out[site.old > 0 & !duplicated(site.old), 1], 
                   out[site.old > 0 & !duplicated(site.old, fromLast = T), 1]), lty = 2)
      abline(v = c(out[out$site > 0 & !duplicated(out$site), 1], 
                   out[out$site > 0 & !duplicated(out$site, fromLast = T), 1]), lwd = 1.5, 
             col = "firebrick")
      axis(1, at = seq(out[1, 1], out[nrow(out), 1], length = 10), 
           labels = format(seq(out[1, 1], out[nrow(out), 1], length = 10), "%d-%b"))
      mtext("Date", 1, outer = T, line = 1.6, cex = 1.2)
      par(opar)
    }
  }
  
  sm        <- as.data.frame(sm)
  names(sm) <- c("site", "Lon", "Lat", "Lon.2.5%", "Lon.40%", "Lon.80%", "Lon.97.25%",  
                 "Lat.2.5%", "Lat.40%", "Lat.80%", "Lat.97.25%")
  
  diff <- c(apply(cbind(out$datetime[-nrow(out)], 
                        out$datetime[-1]), 1, function(x) c(x[2] - x[1])/60/60), 0)
  
  out.gl <- data.frame(tFirst = as.POSIXct("1900-01-01 01:01", "GMT"), tSecond = as.POSIXct("1900-01-01 01:01", "GMT"),  type = 0, site = 0, fixed = 0, diff.max = 0)
  rw <- 1
  for (k in 1:(nrow(out) - 1)) {
    if (as.numeric(difftime(out$datetime[k], out$datetime[k + 1])) < 24 & out$datetime[k] != out$datetime[k + 1]) {
      out.gl[rw, 1] <- out$datetime[k]
      out.gl[rw, 2] <- out$datetime[k + 1]
      out.gl[rw, 3] <- ifelse(out$Rise[k], 1, 2)
      out.gl[rw, 4] <- out$site[k]
      out.gl[rw, 5] <- out$fixed[k]
      out.gl[rw, 6] <- max(diff[k:(k + 1)])
      rw <- rw + 1
    }
  }
  out.gl <- out.gl[out.gl$diff.max < 23,]
  
  list(twl = out.gl[,c(1,2)], site = out.gl$site, summary = sm)
}


##' Filter for unrealistic positions within a track based on distance
##'
##' The filter identifies unrealistic positions. The maximal distance per
##' hour/day can be set corresponding to the particular species.
##'
##' Note that this type of filter significantly depends on the calibration
##' (\code{degElevation}). Especially during equinox periods. In contrast, the
##' (\code{loessFilter}) is independent from positions (uses twilight times) 
##' and therefore superior.
##'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
##' twilight")
##' @param distance the maximal distance in km per \code{units}. Distances above
##' will be considered as unrealistic.
##' @param units the time unite corresponding to the distance. Default is
##' "hour", alternative option is "day".
##' @return Logical \code{vector}. TRUE means the particular position passed the filter.
##' @author Simeon Lisovski, Fraenzi Korner-Nievergelt
##' @importFrom stats loess
##' @export distanceFilter
distanceFilter <- function(tFirst, tSecond, type, twl, degElevation = -6, distance, units = "hour") {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  if(units=="days") units <- distance/24
  
  tFirst <- as.POSIXct(tab$tFirst, tz = "GMT")
  tSecond<- as.POSIXct(tab$tSecond, tz = "GMT")
  type <- tab$type
  tSunTransit <- tFirst + (tSecond-tFirst)/2
  crds  <- coord(tab, degElevation, note=FALSE)
  
  crds[is.na(crds[,2]),2] <- 999
  
  difft <- as.numeric(difftime(tSunTransit[-length(tSunTransit)],tSunTransit[-1],units="hours"))
  diffs <- abs(as.numeric(i.loxodrom.dist(crds[-nrow(crds), 1], crds[-nrow(crds), 2], crds[-1, 1], crds[-1, 2]))/difft)
  
  index <- rep(TRUE,length(tFirst))
  index[diffs>distance] <- FALSE
  index[crds[,2]==999] <- TRUE
  
  
  cat(paste("Note: ",length(index[!index])," of ",length(index[crds[,2]!=999])," positions were filtered (",floor((length(index[!index])*100)/length(index[crds[,2]!=999]))," %)",sep=""))
  index
}


#' Transformation of *.gle files
#'
#' Function to transform *.gle files derived by the software GeoLocator for
#' further analyses in \bold{\code{GeoLight}}.
#'
#' The *.gle file derived by the software "GeoLocator" (Swiss Ornithological
#' Institute) is a table with interpolated light intensities over time gathered
#' from the *.glf file. Furthermore a column defines wether the light intensity
#' passes the defined light intensity threshold in the morning (sunrise) or in
#' the evening (sunset). This information is used in \code{gleTrans()}, to
#' create a table with two subsequent twilight events (\emph{tFirst, tSecond})
#' and \emph{type} defining wether \emph{tFirst} refers to sunrise (1) or
#' sunset (2). Date and time information will be transferred into a
#' \code{GeoLight} appropriate format (see: \code{\link{as.POSIXct}}).
#'
#' @param file the full patch and filename with suffix of the *.gle file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{glfTrans}}
#' @importFrom utils read.table
#' @export gleTrans
gleTrans <- function(file) {
  
  
  gle1 <- read.table(file,sep="\t",skip=16,col.names=c("date","light","1","2","3","4","5","6","7","8","type")) # read file
  gle1 <- subset(gle1,gle1$type>0,select=c("date","type"))
  
  # Date transformation
  year   <- as.numeric(substring(gle1$date,7,10))
  month  <- as.numeric(substring(gle1$date,4,5))
  day    <- as.numeric(substring(gle1$date,1,2))
  hour   <- as.numeric(substring(gle1$date,12,13))
  min    <- as.numeric(substring(gle1$date,15,16))
  gmt.date <- paste(year,"-", month,"-",day," ",hour,":",min,":",0,sep="")
  gmt.date <- as.POSIXct(strptime(gmt.date, "%Y-%m-%d %H:%M:%S"), "UTC")
  
  gle <- data.frame(date=gmt.date,type=gle1$type)
  
  
  opt <- data.frame(tFirst=as.POSIXct("1900-01-01 01:01","UTC"),tSecond=as.POSIXct("1900-01-01 01:01","UTC"),type=0)
  
  row <- 1
  for (i in 1:(length(gmt.date)-1))
  {
    if (abs(as.numeric(difftime(gle$date[i],gle$date[i+1],units="hours")))< 18 & gle$date[i] != gle$date[i+1])
    {
      opt[row,1] <- gle$date[i]
      opt[row,2] <- gle$date[i+1]
      opt[row,3] <- gle$type[i]
      
      row <- row+1
    }
  }
  
  return(opt)
}

#' Transformation of *.glf files
#'
#' Transform *.glf files derived by the software GeoLocator for further
#' analyses in \bold{\code{GeoLight}}.
#'
#' The *.glf files produced by the software GeoLocator (Swiss Ornithological
#' Institute) is a table with light intensity measurements over time.
#' \code{glfTrans} produces a table with these measurements and transfer the
#' data and time information into the format required by \bold{\code{GeoLight}}
#' format (see: \code{\link{as.POSIXct}}).
#'
#' @param file the full patch and filename with suffix of the *.glf file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{gleTrans}}; \code{\link{luxTrans}} for transforming
#' *.lux files produced by \emph{Migrate Technology Ltd}
#' @importFrom utils read.table
#' @export glfTrans
glfTrans <- function(file="/path/file.glf") {
  
  
  glf1 <- read.table(file,sep="\t",skip=6,col.names=c("datetime","light","1","2","3")) # read file
  
  # Date transformation
  year   <- as.numeric(substring(glf1$datetime,7,10))
  month  <- as.numeric(substring(glf1$datetime,4,5))
  day    <- as.numeric(substring(glf1$datetime,1,2))
  hour   <- as.numeric(substring(glf1$datetime,12,13))
  min    <- as.numeric(substring(glf1$datetime,15,16))
  gmt.date <- paste(year,"-", month,"-",day," ",hour,":",min,":",0,sep="")
  gmt.date <- as.POSIXct(strptime(gmt.date, "%Y-%m-%d %H:%M:%S"), "UTC")
  
  glf <- data.frame(datetime=gmt.date,light=glf1$light)
  
  return(glf)
}

#' Hill-Ekstrom calibration
#'
#' Hill-Ekstrom calibration for one or multiple stationary periods.
#'
#' The \emph{Hill-Ekstrom calibration} has been suggested by Hill & Braun
#' (2001) and Ekstrom (2004), and allows for calibrating data during stationary
#' periods at unknown latitudinal positions. The Hill-Ekstrom calibration bases
#' on an increasing error range in latitudes with an increasing mismatch
#' between light level threshold and the used sun angle. This error is strongly
#' amplified with proximity to the equinox times due to decreasing slope of day
#' length variation with latitude. Furthermore, the sign of the error switches
#' at the equinox, i.e. latitude is overestimated before the equinox and
#' underestimated after the equinox (or vice versa depending on autumnal/vernal
#' equinox, hemisphere, and sign of the mismatch between light level threshold
#' and sun angle). When calculating the positions of a stationary period, the
#' variance in latitude is minimal if the sun elevation angle fits to the
#' defined light level threshold. Moreover, the accuracy of positions increases
#' with decreasing variance in latitudes. \bold{However, the method is only
#' applicable for stationary periods and under stable shading intensities}. The
#' plot produced by the function may help to judge visually if the calculated
#' sun elevation angles are realistic (e.g. site 2 in the example below) or not
#' (e.g. site 3 in the example below.
#'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param site a \code{numerical vector} assigning each row to a particular
##' period. Stationary periods in numerical order and values >0,
##' migration/movement periods 0
##' @param start.angle a single sun elevation angle. The combined process of
##' checking for minimal variance in resulting latitude, which is the initial
##' value for the sun elevation angle in the iterative process of identifying
##' the latitudes with the least variance
##' @param distanceFilter logical, if TRUE the \code{\link{distanceFilter}} will
##' be used to filter unrealistic positions
##' @param distance if \code{distanceFilter} is set \code{TRUE} a threshold
##' distance in km has to be set (see: \code{\link{distanceFilter}})
##' @param plot logical, if TRUE the function will give a plot with all relevant
##' information
##' @return A \code{vector} of sun elevation angles corresponding to the
##' Hill-Ekstrom calibration for each defined period.
##' @author Simeon Lisovski
##' @references Ekstrom, P.A. (2004) An advance in geolocation by light.
##' \emph{Memoirs of the National Institute of Polar Research}, Special Issue,
##' \bold{58}, 210-226.
##'
##' Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step:
##' Latitude. In: \emph{Electronic Tagging and Tracking in Marine Fisheries}
##' (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The
##' Netherlands.
##'
##' Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F.,
##' Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
##' precision affected by environmental factors. \emph{Methods in Ecology and
##' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
##' @examples
##'
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' residency <- with(hoopoe2, changeLight(tFirst,tSecond,type, rise.prob=0.1, 
##'                   set.prob=0.1, plot=FALSE, summary=FALSE))
##' HillEkstromCalib(hoopoe2,site = residency$site)
##'
##' @importFrom grDevices grey.colors
##' @importFrom graphics abline axis mtext legend lines par plot points text
##' @importFrom stats var na.omit
##' @export HillEkstromCalib
HillEkstromCalib <- function(tFirst, tSecond, type, twl, site, start.angle=-6, distanceFilter=FALSE, distance, plot=TRUE) {
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  tFirst <- tab$tFirst
  tSecond <- tab$tSecond
  type <- tab$type
  
  sites <- as.numeric(length(levels(as.factor(site[as.numeric(site)!=0]))))
  
  HECalib <- rep(NA,sites)
  
  
  for(j in 1:sites) {
    b <- 0
    start <- start.angle
    
    repeat{
      # backwards
      if(start-((b*0.1)-0.1) < -9) {
        HECalib[j] <- NA
        break
      }
      t0 <- var(na.omit(coord(tab[site==j,], degElevation = start-((b*0.1)-0.1), note = F)[,2]))
      t1 <- var(na.omit(coord(tab[site==j,], degElevation = start-(b*0.1), note = F)[,2]))
      t2 <- var(na.omit(coord(tab[site==j,], degElevation = start-((b*0.1)+0.1), note = F)[,2]))
      if(sum(is.na(c(t0,t1,t2)))>0) {
        HECalib[j] <- NA
        break
      }
      if(t0>t1 & t1<t2) {
        HECalib[j] <- start-(b*0.1)
        break}
      
      # forward
      if(start-((b*0.1)-0.1) > 9) {
        HECalib[j] <- NA
        break
      }
      f0 <- var(na.omit(coord(tab[site==j,], degElevation = start+((b*0.1)-0.1), note = F)[,2]))
      f1 <- var(na.omit(coord(tab[site==j,], degElevation = start+(b*0.1), note = F)[,2]))
      f2 <- var(na.omit(coord(tab[site==j,], degElevation = start+((b*0.1)+0.1), note = F)[,2]))
      if(sum(is.na(c(f0,f1,f2)))>0) {
        HECalib[j] <- NA
        break
      }
      if(f0>f1 & f1<f2) {
        HECalib[j] <- start+(b*0.1)
        break}
      
      b <- b+1
    }
  }
  
  if(plot) {
    
    opar <- par(mfrow = c(sites, 2), oma = c(3, 5, 6, 2), mar = c(2, 2, 2, 2))
    for(j in 1:sites){
      
      if(is.na(HECalib[j])){
        
        plot(0,0,cex=0,pch=20,col="white",ylab="",xlab="",xaxt="n",yaxt="n", bty = "n")
        text(0,0,"NA",cex=2)
        plot(0,0,cex=0,pch=20,col="white",ylab="",xlab="",xaxt="n",yaxt="n", bty="n")
        
      } else {
        
        angles <- c(seq(HECalib[j]-2,HECalib[j]+2,0.2))
        
        latM <- matrix(ncol=length(angles),nrow=length(tFirst[site==j]))
        
        for(i in 1:ncol(latM)){
          latM[,i] <- coord(tab[site==j,], degElevation = c(angles[i]),note=F)[,2]
        }
        
        latT <- latM
        var1 <- rep(NA,ncol(latT))
        n1   <- rep(NA,ncol(latT))
        min  <- rep(NA,ncol(latT))
        max  <- rep(NA,ncol(latT))
        for(t in 1:length(var1)){
          var1[t] <- var(na.omit(latT[,t]))
          n1[t]   <- length(na.omit(latT[,t]))
          min[t]  <- if(length(na.omit(latT[,t]))<=1) NA else min(na.omit(latT[,t]))
          max[t]  <- if(length(na.omit(latT[,t]))<=1) NA else max(na.omit(latT[,t]))
        }
        
        
        colors <- grey.colors(length(angles))
        plot(tFirst[site==j],latT[,1],ylim=c(min(na.omit(min)),max(na.omit(max))),type="o",cex=0.7,pch=20,col=colors[1],ylab="",xlab="")
        for(p in 2:ncol(latT)){
          lines(tFirst[site==j],latM[,p],type="o",cex=0.7,pch=20,col=colors[p])
        }
        lines(tFirst[site==j],coord(tab[site==j,], degElevation = HECalib[j], note=F)[,2],col="tomato2",type="o",lwd=2,cex=1,pch=19)
        if(j==sites) mtext("Latitude",side=2,line=3)
        if(j==sites) mtext("Date",side=1,line=2.8)
        
        plot(angles,var1,type="o",cex=1,pch=20,ylab="")
        lines(angles,var1,type="p",cex=0.5,pch=7)
        abline(v=HECalib[j],lty=2,col="red",lwd=1.5)
        par(new=T)
        plot(angles,n1,type="o",xaxt="n",yaxt="n",col="blue",pch=20,ylab="")
        points(angles,n1,type="p",cex=0.5,pch=8,col="blue")
        axis(4)
        if(j==sites) mtext("Sun elevation angles",side=1,line=2.8)
        legend("topright",c(paste(HECalib[j]," degrees",sep=""),"sample size","variance"),pch=c(-1,20,20),lty=c(2,2,2),lwd=c(1.5,0,0),col=c("red","blue","black"),bg="White")
      }
    }
    
    mtext("Hill-Ekstrom Calibration",cex=1.5,line=0.8,outer=T)
    par(opar)
    
  }
  
  
  return(HECalib)
}



i.JC2000 <- function(jD) {
  
  #--------------------------------------------------------------------------------------------------------
  # jD: julian Date
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  
  jC<- (jD - 2451545)/36525
  
  return(jC)
  
}

i.deg <- function(Rad) {
  
  #------------------------------------------------------------
  # Deg: 	The input angle in radian
  #------------------------------------------------------------
  
  options(digits=10)
  
  
  return(Rad * (180/pi))
  
}

i.frac <- function(In) {
  
  #------------------------------------------------------------
  # In: 	numerical Number
  #------------------------------------------------------------
  
  options(digits=10)
  
  
  return(In - floor(In))
  
}

i.get.outliers<-function(residuals, k=3) {
  x <- residuals
  # x is a vector of residuals
  # k is a measure of how many interquartile ranges to take before saying that point is an outlier
  # it looks like 3 is a good preset for k
  QR<-quantile(x, probs = c(0.25, 0.75))
  IQR<-QR[2]-QR[1]
  Lower.band<-QR[1]-(k*IQR)
  Upper.Band<-QR[2]+(k*IQR)
  delete<-which(x<Lower.band |  x>Upper.Band)
  return(as.vector(delete))
}

i.julianDate <- function(year,month,day,hour,min) {
  
  #--------------------------------------------------------------------------------------------------------
  # Year:		Year as numeric e.g. 2010
  # Month:	Month as numeric e.g.   1
  # Day:		Day as numeric e.g.		1
  # Hour:		Hour as numeric e.g.   12
  # Min:		Minunte as numeric e.g. 0
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=15)
  
  fracOfDay	<- hour/24 + min/1440
  
  # Julian date (JD)
  # ------------------------------------
  
  index1 <- month <= 2
  if(sum(index1) > 0)
  {
    year[index1]  <- year[index1] -1
    month[index1] <- month[index1] +12
  }
  
  index2 <- (year*10000)+(month*100)+day <= 15821004
  
  JD <- numeric(length(year))
  if(sum(index2)>0)
  {
    JD[index2] <- floor(365.25*(year[index2]+4716)) + floor(30.6001*(month[index2]+1)) + day[index2] + fracOfDay[index2] - 1524.5
  }
  index3 <- year*10000+month*100+day >= 15821015
  if (sum(index3)>0)
  {
    a <- floor(year/100)
    b <- 2 - a + floor(a/4)
    
    JD[index3] <- floor(365.25*(year[index3]+4716)) + floor(30.6001*(month[index3]+1)) + day[index3] + fracOfDay[index3] + b[index3] - 1524.5
  }
  JD[!index2&!index3] <- 1
  
  return(JD)
  
}

i.loxodrom.dist <- function(x1, y1, x2, y2, epsilon=0.0001){
  dis<-numeric(length(x1))
  rerde<-6368
  deltax<-abs(x2*pi/180-x1*pi/180)
  deltay<-abs(y2*pi/180-y1*pi/180)
  tga<-deltax/(log(tan(pi/4+y2*pi/360))-log(tan(pi/4+y1*pi/360)))
  
  dis[abs(x1-x2)<epsilon&abs(y1-y2)<epsilon]<-0
  dis[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]<-abs(cos(y1[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)]*pi/180)*deltax[abs(y1-y2)<epsilon&(abs(x1-x2)>epsilon)])
  dis[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]<-abs(deltay[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]/cos((pi-atan(tga[(tga<0)&(abs(x1-x2)>epsilon)&(abs(y1-y2)>epsilon)]))))
  dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]<--deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y1-y2>epsilon)]))
  dis[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]<-abs(deltay[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(tga>=0)&(abs(x1-x2)>epsilon)&(y2-y1>epsilon)])))
  dis[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y2-y1>epsilon)])))
  dis[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]<-abs(deltay[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)]/cos(atan(tga[(abs(x1-x2)<epsilon)&(y1-y2>epsilon)])))
  dis*rerde
}

i.preSelection <- function(datetime, light, LightThreshold){
  
  dt <- cut(datetime,"1 hour")
  st <- as.POSIXct(levels(dt),"UTC")
  
  raw <- data.frame(datetime=dt,light=light)
  
  h  <- tapply(light,dt,max)
  df1 <- data.frame(datetime=st+(30*60),light=as.numeric(h))
  
  smooth <- i.twilightEvents(df1[,1], df1[,2], LightThreshold)
  smooth <- data.frame(id=1:nrow(smooth),smooth)
  raw    <- i.twilightEvents(datetime, light, LightThreshold)
  raw <- data.frame(id=1:nrow(raw),raw)
  
  ind2 <- rep(NA,nrow(smooth))
  for(i in 1:nrow(smooth)){
    tmp <- subset(raw,datetime>=(smooth[i,2]-(90*60)) & datetime<=(smooth[i,2]+(90*60)))
    
    if(smooth[i,3]==1) ind3 <- tmp$id[which.min(tmp[,2])]
    if(smooth[i,3]==2) ind3 <- tmp$id[which.max(tmp[,2])]
    ind2[i] <- ind3
  }
  
  
  res <- data.frame(raw,mod=1)
  res$mod[ind2] <- 0
  
  return(res)
}

i.rad <- function(Deg) {
  
  #------------------------------------------------------------
  # Deg: 	The input angle in degrees
  #------------------------------------------------------------
  
  options(digits=10)
  
  
  return(Deg * (pi/180))
  
}

i.radDeclination <- function(radEclipticalLength,radObliquity) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane.
  # EadObliquity:
  #-------------------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  dec <- asin(sin(radObliquity)*sin(radEclipticalLength))
  
  return(dec)
  
}

i.radEclipticLongitude <- function(jC) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00
  #-------------------------------------------------------------------------------------------------------------------
  
  
  
  options(digits=10)
  
  radMeanAnomaly <- 2*pi*i.frac(0.993133 + 99.997361*jC)
  EclipticLon    <- 2*pi*i.frac(0.7859452 + radMeanAnomaly/(2*pi) + (6893*sin(radMeanAnomaly) + 72*sin(2*radMeanAnomaly) + 6191.2*jC) / 1296000)
  
  return(EclipticLon)
  
}

i.radGMST <- function(jD,jD0,jC,jC0) {
  
  #--------------------------------------------------------------------------------------------------------
  # jD:  Julian Date with Hour and Minute
  # jD0: Julan Date at t 0 UT
  # jC:  Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
  # jC0: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00) at t 0 UT
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  
  UT  <- 86400 * (jD-jD0)
  st0 <- 24110.54841 + 8640184.812866*jC0 + 1.0027379093*UT + (0.093104 - 0.0000062*jC0)*jC0*jC0
  gmst<- (((2*pi)/86400)*(st0%%86400))
  
  return(gmst)
  
}

i.radObliquity <- function(jC) {
  
  #--------------------------------------------------------------------------------------------------------
  # jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
  #--------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  
  degObliquity <- 23.43929111 - (46.8150 + (0.00059 - 0.001813 *jC)*jC) *jC/3600
  radObliquity <- i.rad(degObliquity)
  
  return(radObliquity)
  
}

i.radRightAscension <- function(RadEclipticalLength,RadObliquity) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane.
  # EadObliquity:
  #-------------------------------------------------------------------------------------------------------------------
  
  options(digits=10)
  
  index1 <- (cos(RadEclipticalLength) < 0)
  
  res <- numeric(length(RadEclipticalLength))
  if (sum(index1)>0)
  {
    res[index1] <- (atan((cos(RadObliquity[index1])*sin(RadEclipticalLength[index1]))/cos(RadEclipticalLength[index1])) + pi)
  }
  
  index2 <- (cos(RadEclipticalLength) >= 0)
  if (sum(index2)>0)
  {
    res[index2] <- (atan((cos(RadObliquity[index2])*sin(RadEclipticalLength[index2]))/cos(RadEclipticalLength[index2])))
  }
  
  
  return(res)
  
}

i.setToRange <- function(Start,Stop,Angle) {
  
  #-------------------------------------------------------------------------------------------------------------------
  # Start:	 Minimal value of the range in degrees
  # Stop: 	 Maximal value of the range in degrees
  # Angle:	 The angle that should be fit to the range
  #-------------------------------------------------------------------------------------------------------------------
  
  options(digits=15)
  
  angle <- Angle
  range <- Stop - Start
  
  
  
  index1 <- angle >= Stop
  if (sum(index1, na.rm = T)>0) angle[index1] <- angle[index1] - (floor((angle[index1]-Stop)/range)+1)*range
  
  index2 <- angle < Start
  if(sum(index2, na.rm = T)>0) angle[index2] <- angle[index2]  + ceiling(abs(angle[index2] -Start)/range)*range
  
  return(angle)
  
}

i.sum.Cl <- function(object) {
  
  if(all(names(object)%in%c("riseProb","setProb","rise.prob","set.prob","site","migTable"))){
    cat("\n")
    cat("Probability threshold(s):")
    cat(rep("\n",2))
    if(!is.na(object$rise.prob)) cat(paste("	Sunrise: ",object$rise.prob))
    if(!is.na(object$set.prob)) cat(paste("	Sunset: ",object$set.prob))
    cat(rep("\n",3))
    cat("Migration schedule table:")
    cat(rep("\n",2))
    
    print(object$migTable,quote=FALSE)
  } else {
    cat("Error: List must be the output list of the changeLight function.")
  }
}


gl.loglik <- function(crds, Twilight, Rise, degElevation, parms, method, twilight = NULL) {
  t.tw <- twilight(Twilight, lon = crds[1], lat = crds[2], 
                   rise = Rise, zenith = 90 - degElevation, 
                   iters = 6)
  
  diff.sr <- as.numeric(difftime(Twilight[Rise], t.tw[Rise], units = "mins"))
  diff.ss <- as.numeric(difftime(t.tw[!Rise], Twilight[!Rise], units = "mins"))
  
  if(method=="log-norm") {
    if(is.null(twilight)) {
      ll <- sum(dlnorm(c(diff.sr, diff.ss), parms[1], parms[2], log = F), na.rm = T)
    } else {
      if(twilight=="sr"){
        ll <- sum(dlnorm(diff.sr, parms[1], parms[2], log = F), na.rm = T)
      }
      if(twilight=="ss") {
        ll <- sum(dlnorm(diff.ss, parms[1], parms[2], log = F), na.rm = T)
      }}
  } else {
    if(is.null(twilight)) {
      ll <- sum(dgamma(c(diff.sr, diff.ss), parms[1], parms[2], log = F), na.rm = T)
    } else {
      if(twilight=="sr"){
        ll <- sum(dgamma(diff.sr, parms[1], parms[2], log = F), na.rm = T)
      }
      if(twilight=="ss") {
        ll <- sum(dgamma(diff.ss, parms[1], parms[2], log = F), na.rm = T)
      }}
  }
  
  return(ifelse(!is.infinite(abs(ll)), ll, -10000))
}


#' Filter to remove noise in light intensity measurements during the night
#'
#' The filter identifies and removes light intensities oczillating around the
#' baseline or few light intensities resulting in a short light peak during the
#' night. Such noise during the night will increase the calculated twilight
#' events using the function \code{\link{twilightCalc}} and therewith the
#' manual work to remove these false twilight events.
#'
#' The filter searches for light levels above the baseline and compares the
#' prior and posterior levels. If these values are below the threshold the
#' particular light level will be reduced to the baseline. A few (usually two)
#' iterations might be enough to remove most noise during the night (however,
#' not if such noise occurs at the begining or at the end were not enough prior
#' or posterior values are available).
#'
#' @param light \code{numerical} value of the light intensity (usually
#' arbitrary units).
#' @param baseline the light intensity baseline (no light). If \code{Default},
#' it will be calculated as the most frequent value below the mean light
#' intensities.
#' @param iter a \code{numerical} value, specifying how many iterations should
#' be computed (see details).
#' @return numerical \code{vector} with the new light levels. Same length as
#' the initial light vector.
#' @author Simeon Lisovski
#' @examples
#'
#' night <- rep(0,50); night[runif(4,0,50)] <- 10; night[runif(4,0,50)] <- -5
#' nightday <- c(night,rep(30,50))
#' plot(nightday,type="l",ylim=c(-5,30),ylab="light level",xlab="time (time)")
#' light2 <- lightFilter(nightday, baseline=0, iter=4)
#' lines(light2,col="red")
#' legend("bottomright",c("before","after"),lty=c(1,1),col=c("black","red"),bty="n")
#'
#' @export lightFilter
lightFilter <- function(light, baseline=NULL, iter=2){
  
  r <- as.data.frame(table(light))
  r[,1] <- as.character(r[,1])
  nr <- as.numeric(which.max(r$Freq[as.numeric(r[,1])<mean(light)]))
  LightThreshold <- ifelse(is.null(baseline),as.numeric(r[nr,1]),baseline)
  
  light[light<LightThreshold] <- LightThreshold
  
  
  index <- which(light<mean(light) & light!= LightThreshold)
  
  rep   <- rep(FALSE,length(light))
  
  for(i in 1:iter){
    
    for(i in index[index>5 & index<(length(light)-5)]){
      
      back=FALSE
      if(any(light[seq(i-5,i)]==LightThreshold)) back <- TRUE
      forw=FALSE
      if(any(light[seq(i,i+5)]==LightThreshold)) forw <- TRUE
      
      if(back & forw) rep[i] <- TRUE
      
    }
    
    light[rep] <- LightThreshold
    
  }
  
  light
}


##' Filter to remove outliers in defined twilight times based on smoother function
##'
##' This filter defines outliers based on residuals from a local polynomial
##' regression fitting provcess (\code{\link{loess}}).
##'
##'
##' @param tFirst vector of sunrise/sunset times (e.g. 2008-12-01 08:30).
##' @param tSecond vector of of sunrise/sunset times (e.g. 2008-12-01 17:30).
##' @param type vector of either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively.
##' @param twl data.frame containing twilights and at least \code{tFirst}, \code{tSecond} and \code{type} (alternatively give each parameter separately).
##' @param k a measure of how many interquartile ranges to take before saying
##' that a particular twilight event is an outlier
##' @param plot codelogical, if TRUE a plot indicating the filtered times will
##' be produced.
##' @return Logical \code{vector} matching positions that pass the filter.
##' @author Simeon Lisovski & Eldar Rakhimberdiev
##' @importFrom stats loess predict residuals
##' @importFrom ggplot2 ggplot geom_point labs theme_bw theme element_blank
##' @importFrom patchwork wrap_plots
##' @export loessFilter

loessFilter <- function(tFirst, tSecond, type, twl, k = 3, plot = TRUE){
  
  tab <- i.argCheck(as.list(environment())[sapply(environment(), FUN = function(x) any(class(x)!='name'))])   
  
  tw <- data.frame(datetime = .POSIXct(c(tab$tFirst, tab$tSecond), "GMT"), 
                   type = c(tab$type, ifelse(tab$type == 1, 2, 1)))
  tw <- tw[!duplicated(tw$datetime),]
  tw <- tw[order(tw[,1]),]
  
  hours <- as.numeric(format(tw[,1],"%H"))+as.numeric(format(tw[,1],"%M"))/60
  
  for(t in 1:2){
    cor <- rep(NA, 24)
    for(i in 0:23){
      cor[i+1] <- max(abs((c(hours[tw$type==t][1],hours[tw$type==t])+i)%%24 -
                            (c(hours[tw$type==t],hours[tw$type==t][length(hours)])+i)%%24),na.rm=T)
    }
    hours[tw$type==t] <- (hours[tw$type==t] + (which.min(round(cor,2)))-1)%%24
  }
  
  dawn <- data.frame(id=1:sum(tw$type==1),
                     datetime=tw$datetime[tw$type==1],
                     type=tw$type[tw$type==1],
                     hours = hours[tw$type==1], filter=FALSE)
  
  dusk <- data.frame(id=1:sum(tw$type==2),
                     datetime=tw$datetime[tw$type==2],
                     type=tw$type[tw$type==2],
                     hours = hours[tw$type==2], filter=FALSE)
  
  
  for(d in seq(30,k,length=5)){
    
    predict.dawn <- predict(loess(dawn$hours[!dawn$filter]~as.numeric(dawn$datetime[!dawn$filter]),span=0.1))
    predict.dusk <- predict(loess(dusk$hours[!dusk$filter]~as.numeric(dusk$datetime[!dusk$filter]),span=0.1))
    
    del.dawn <-	i.get.outliers(as.vector(residuals(loess(dawn$hours[!dawn$filter]~
                                                           as.numeric(dawn$datetime[!dawn$filter]),span=0.1))),k=d)
    del.dusk <-	i.get.outliers(as.vector(residuals(loess(dusk$hours[!dusk$filter]~
                                                           as.numeric(dusk$datetime[!dusk$filter]),span=0.1))),k=d)
    
    if(length(del.dawn)>0) dawn$filter[!dawn$filter][del.dawn] <- TRUE
    if(length(del.dusk)>0) dusk$filter[!dusk$filter][del.dusk] <- TRUE
  }
  
  if(plot){
    
    sunrise_plot <- ggplot() +
      geom_point(aes(dawn$datetime[dawn$type==1],
                     dawn$hours[dawn$type==1]), 
                 shape = 3) + 
      geom_line(aes(dawn$datetime[!dawn$filter],
                    predict(loess(dawn$hours[!dawn$filter]~as.numeric(dawn$datetime[!dawn$filter]),
                                  span=0.1)))) + 
      geom_point(aes(dawn$datetime[dawn$filter],
                     dawn$hours[dawn$filter]), 
                 col = "red", 
                 shape = 3 ) + 
      labs(title = "Sunset/sunrise hours (rescaled)",
           y = "Sunrise",
           x = "") +
      theme_bw() +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank() ,
            axis.text.x = element_blank()) 
    
    sunset_plot <- ggplot() +
      geom_point(aes(dusk$datetime[dusk$type==2], 
                     dusk$hours[dusk$type==2]), 
                 shape = 3) + 
      geom_line(aes(dusk$datetime[!dusk$filter],
                    predict(loess(dusk$hours[!dusk$filter]~as.numeric(dusk$datetime[!dusk$filter]),
                                  span=0.1)))) + 
      geom_point(aes(dusk$datetime[dusk$filter],
                     dusk$hours[dusk$filter]), 
                 col = "red", 
                 shape = 3) + 
      labs(x = "Time",
           y = "Sunset") +
      theme_bw() +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank())
    
    print(wrap_plots(sunrise_plot / sunset_plot))
    
  } 
  all <- rbind(subset(dusk,filter),subset(dawn,filter))
  
  filter <- rep(FALSE,length(tab$tFirst))
  filter[as.POSIXct(tab$tFirst, "GMT")%in%all$datetime | as.POSIXct(tab$tSecond, "GMT")%in%all$datetime] <- TRUE
  
  return(!filter)
}


#' Transformation of *.lux files
#'
#' Transform *.lux files derived from \emph{Migrate Technology Ltd} geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' The *.lux files produced by \emph{Migrate Technology Ltd} are table with
#' light intensity measurements over time. \code{luxTrans} produces a table
#' with these measurements and transfer the data and time information into the
#' format required by \bold{\code{GeoLight}} format (see:
#' \code{\link{as.POSIXct}}).
#'
#' @param file the full patch and filename with suffix of the *.lux file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{gleTrans}} for transforming *.glf files produced by the
#' software GeoLocator (\emph{Swiss Ornithological Institute})
#' @importFrom utils read.table
#' @export luxTrans
luxTrans <- function(file) {
  
  lux1 <- read.table(file,sep="\t",skip=21,col.names=c("datetime","time")) # read file
  lux <- data.frame(datetime=as.POSIXct(strptime(lux1[,1],format="%d/%m/%Y %H:%M:%S",tz="UTC")),light=lux1[,2])
  
  return(lux)
}


#' Transformation of *.lig files
#'
#' Transform *.lig files derived from geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' @param file the full patch and filename with suffix of the *.lux file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Tamara Emmenegger
#' @seealso \code{\link{gleTrans}} for transforming *.glf files produced by the
#' software GeoLocator (\emph{Swiss Ornithological Institute})
#' @importFrom utils read.csv
#' @export ligTrans
ligTrans <- function(file) {
  df <- cbind(read.csv(file,row.names=NULL)[,c(2,4)])
  colnames(df) <- c("datetime", "light")
  df$datetime <- as.POSIXct(strptime(df$datetime,format="%d/%m/%y %H:%M:%S",tz="UTC"))
  return(df)
}

#' Transformation of staroddi files
#'
#' Transform staroddi files derived from geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' @param file the full patch and filename of the staroddi file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Tamara Emmenegger
#' @importFrom utils read.table
#' @export staroddiTrans
staroddiTrans <- function(file){
  df <- cbind(read.table(file,sep="\t")[,c(2,4)])
  colnames(df) <- c("datetime", "light")
  df$datetime <- as.POSIXct(strptime(df$datetime,format="%d/%m/%y %H:%M:%S",tz="UTC"))
  return(df)
}


#' Transformation of *.trn files
#'
#' Transform *.trn files derived from geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#'
#' @param file the full patch and filename of the *.trn file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Tamara Emmenegger
#' @importFrom utils read.table
#' @export trnTrans
trnTrans<-function(file){
  data<-read.table(file,sep=",")
  tFirst <- vector("numeric",length=(length(data[,1])-1))
  tSecond <- vector("numeric",length=(length(data[,1])-1))
  type <- vector("numeric",length=(length(data[,1])-1))
  for (i in 1:(length(data[,1])-1)){
    date1<-as.Date(substr(as.character(data$V1[i]),1,8),format="%d/%m/%y")
    tFirst[i] <- as.POSIXct(paste(as.character(date1),prefix=substr(as.character(data$V1[i]),10,17)),tz="UTC")
    date2<-as.Date(substr(as.character(data$V1[i+1]),1,8),format="%d/%m/%y")
    tSecond[i] <- as.POSIXct(paste(as.character(date2),prefix=substr(as.character(data$V1[i+1]),10,17)),tz="UTC")
    if(as.character(data$V2[i])=="Sunrise") type[i] <- 1
    if(as.character(data$V2[i])=="Sunset") type[i] <- 2
  }
  output <- data.frame(tFirst=as.POSIXlt(tFirst,origin="1970-01-01",tz="UTC"),tSecond=as.POSIXlt(tSecond,origin="1970-01-01",tz="UTC"),type=type)
  return(output)
}



## ' Summary of migration/movement pattern
##'
##' Function for making a data frame summarising residency and movement pattern.
##'
##' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
##' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
##' @param site a \code{vector}, indicating the residency period of a particular
##' day (see output: \code{\link{changeLight}})
##' @return A \code{data.frame} with end and start date (yyyy-mm-dd hh:mm, UTC)
##' for each stationary period.
##' @author Simeon Lisovski
##' @examples
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' residency <- changeLight(hoopoe2, rise.prob=0.1, set.prob=0.1, plot=FALSE, summary=FALSE)
##' schedule(hoopoe2[,1], hoopoe2[,2], site = residency$site)
##' @export schedule
schedule <- function(tFirst, tSecond, site) {    
  tm <- tFirst + (tSecond - tFirst)/2
  
  site <- ifelse(is.na(site), 0, site)
  
  arr <- tm[which(!is.na(site) & !duplicated(site) & site>0)]
  dep <- tm[which(!is.na(site) & !duplicated(site, fromLast = T) & site>0)]
  
  out <- data.frame(Site =  letters[1:length(arr)], Arrival = arr, Departure = dep)
  if(!is.na(site[1])) out$Arrival[1] <- NA
  if(!is.na(site[length(tFirst)])) out$Departure[nrow(out)] <- NA
  
  out
}


i.twilightEvents <- function (datetime, light, LightThreshold) 
{
  df <- data.frame(datetime, light)
  ind1 <- which((df$light[-nrow(df)] < LightThreshold & df$light[-1] > 
                   LightThreshold) | (df$light[-nrow(df)] > LightThreshold & 
                                        df$light[-1] < LightThreshold) | df$light[-nrow(df)] == 
                  LightThreshold)
  bas1 <- cbind(df[ind1, ], df[ind1 + 1, ])
  bas1 <- bas1[bas1[, 2] != bas1[, 4], ]
  x1 <- as.numeric(unclass(bas1[, 1]))
  x2 <- as.numeric(unclass(bas1[, 3]))
  y1 <- bas1[, 2]
  y2 <- bas1[, 4]
  m <- (y2 - y1)/(x2 - x1)
  b <- y2 - (m * x2)
  xnew <- (LightThreshold - b)/m
  type <- ifelse(bas1[, 2] < bas1[, 4], 1, 2)
  res <- data.frame(datetime = as.POSIXct(xnew, origin = "1970-01-01", 
                                          tz = "UTC"), type)
  return(res)
}

#' Draws sites of residency and adds a convex hull
#'
#' Draw a map (from the \code{R} Package \code{maps}) showing the defined
#' stationary sites
#'
#'
#' @param crds a \code{SpatialPoints} or \code{matrix} object, containing x
#' and y coordinates (in that order).
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0.
#' @param type either \emph{points}, or \emph{cross} to show all points for each site or only show the mean position of the site with standard deviation.
#' @param quantiles the quantile of the error bars (\emph{cross}) around the median.
#' @param xlim Longitude limits of map.
#' @param ylim Lattitude limits of map.
#' @param hull \code{logical}, if TRUE a convex hull will be plotted around the points of each site.
#' @param palette Color palette for sites.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @details Standard graphical paramters like \code{pch}, \code{cex}, \code{lwd}, \code{lty} and \code{col} are implemented. 
#' The color can be specified as either a vector of colors (e.g. c("blue", "red", ...)) or as a character string indicating a color ramp (at the moment only "random" and "rainbow" is available )
#' @author Simeon Lisovski & Tamara Emmenegger
#' @examples
#' \donttest{
#' data(hoopoe2)
#' hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
#' hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
#' crds <- coord(hoopoe2, degElevation = -6)
#' filter <- distanceFilter(hoopoe2, distance = 30)
#' site <- changeLight(hoopoe2, rise.prob = 0.1, set.prob = 0.1, plot = FALSE, 
#' summary = FALSE)$site
#'siteMap(crds[filter,], site[filter], linewidth=2, shape=20, size=0.5, main="hoopoe2")
#'}
#' @importFrom grDevices chull col2rgb rainbow rgb
#' @importFrom rnaturalearth ne_countries
#' @importFrom ggplot2 theme_bw coord_cartesian labs geom_sf coord_sf geom_point geom_segment geom_polygon guides guide_legend scale_color_hue scale_x_continuous scale_color_manual 
#' @importFrom dplyr %>% mutate as_tibble filter pull group_split .data
#' @importFrom sf sf_use_s2 st_intersection st_bbox st_as_sfc st_transform st_cast st_convex_hull st_union
#' @export siteMap
siteMap <- function(crds, site, type = "points", quantiles = c(0.25, 0.75), xlim = NULL, ylim = NULL, hull = TRUE, palette = 'rainbow', ...) {  
  
  args <- list(...)
  
  if(nrow(crds)!=length(site)) stop("The number of coordinates does not match length of site vector.")
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  dat   <- crds %>% as_tibble() %>% mutate(site=site) %>%
    filter(!is.nan(.data$lat) & site>0)
  
  if(is.null(xlim) | is.null(ylim)) {
    word_sub <- world %>% 
      st_make_valid() %>% 
      st_intersection(st_bbox(dat %>% 
                                st_as_sf(coords = c('lon', 'lat'), crs = 4326)) %>%
                        st_as_sfc())  %>% 
      suppressWarnings() %>% 
      suppressMessages()
  } else {
    word_sub <- world %>% 
      st_make_valid() %>% 
      st_intersection(st_bbox(c(xmin = xlim[1], xmax = xlim[2],
                                ymin = ylim[1], ymax = ylim[2]), crs = 4326) %>%
                        st_as_sfc()) %>% 
      suppressWarnings() %>% 
      suppressMessages()
  }
  
  
  colorSites <- do.call(palette, args = list(n = length(unique(dat$site))))
  
  #### ToDo: Theme - axis not more than map
  base <- ggplot() +
    theme_bw() +
    labs(x = ifelse(sum(names(args) %in% "xlab") == 1, args$xlab, "Longitude"),
         y = ifelse(sum(names(args) %in% "ylab") == 1, args$ylab, "Latitude"),
         title = ifelse(sum(names(args) %in% "main") == 1, args$main, "")) +
    geom_sf(data = word_sub, fill = "lightgray", color = "white") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0))
  
  if(type == "points") {
    
    dat_sf <- dat %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326)
    
    gg <- base +
      geom_sf(data = dat_sf, mapping = aes(geometry = .data$geometry, color = as.factor(site)),
              size = ifelse(any(names(args) == "size"),
                            args$size, 0.5),
              shape = ifelse(any(names(args) == "shape"), args$shape, 16),
              show.legend = ifelse(any(names(args) == "show.legend"), args$show.legend, FALSE)) +
      scale_color_manual(values = colorSites)
  }
  
  if(type == "cross") {
    
    dat_cross <- dat %>% group_split(site) %>%
      lapply(., function(x) {
        x[,1:2] %>% apply(., 2, function(y) quantile(y, probs =  c(quantiles, 0.5), na.rm = T)) %>% 
          as_tibble() %>% mutate(probs =  c(quantiles, 0.5), site = (x %>% pull(site))[1])
      }) %>% Reduce('rbind', .)
    
    dat_pts <- dat_cross %>% filter(.data$probs == 0.5) %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326)
    dat_lns <- dat_cross %>% group_split(site) %>% 
      lapply(., function(x) {
        c(c(st_point(as.numeric(c(x[1,1], x[3,2]))),
            st_point(as.numeric(c(x[2,1], x[3,2])))) %>% st_cast('LINESTRING'),
          c(st_point(as.numeric(c(x[3,1], x[1,2]))),
            st_point(as.numeric(c(x[3,1], x[2,2])))) %>% st_cast('LINESTRING')) %>%
          st_sfc() %>% st_sf() %>% mutate(site = x$site[1])
      }) %>% Reduce('rbind',.) %>% st_set_crs(4326)
    
    gg <- base +
      geom_sf(data = dat_pts, mapping = aes(geometry = .data$geometry, color = as.factor(site)),
              size = ifelse(any(names(args) == "size"),
                            args$size, 0.5),
              shape = ifelse(any(names(args) == "shape"), args$shape, 16),
              show.legend = ifelse(any(names(args) == "show.legend"), args$show.legend, FALSE)) +
      geom_sf(data = dat_lns, mapping = aes(geometry = .data$geometry, color = as.factor(site)),
              linewidth = ifelse(any(names(args) == "linewidth"),
                                 args$linewidth, 0.5),
              linetype = ifelse(any(names(args) == "linetype"), args$linetype, 1),
              show.legend = ifelse(any(names(args) == "show.legend"), args$show.legend, FALSE)) +
      scale_color_manual(values = colorSites)
    
  }
  
  if(hull) {
    
    dat_sf <- dat %>% st_as_sf(coords = c('lon', 'lat'), crs = 4326)
    
    dat_hull <- dat_sf %>%
      group_split(site) %>%
      lapply(., function(x) {
        st_convex_hull(st_union(x[,2])) %>% 
          st_sfc() %>% 
          st_sf() %>% 
          mutate(site = x$site[1])
      })  %>% 
      Reduce('rbind',.) %>% suppressMessages()
    #st_set_crs(4326)
    
    gg <- gg +
      geom_sf(data = dat_hull, 
              mapping = aes(geometry = .data$geometry, color = as.factor(site)),
              alpha = ifelse(any(names(args) == "alpha"), args$alpha, 0.1), 
              linewidth = ifelse(any(names(args) == "linewidth"), args$linewidth, 0.5),
              linetype = ifelse(any(names(args) == "linetype"), args$linetype, 1),
              show.legend = ifelse(any(names(args) == "show.legend"), args$show.legend, FALSE)) 
    
  }
  
  print(gg)
}

##' Write a file which plots a trip in Google Earth
##'
##' This function creates a .kml file from light intensity measurements over
##' time that can ve viewed as a trip in Google Earth.
##'
##'
##' @param file A character expression giving the whole path and the name of the
##' resulting output file including the .kml extension.
##' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
##' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
##' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
##' respectively
##' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
##' twilight"). Either a single value, a \code{vector} with the same length as
##' \code{tFirst}.
##' @param col.scheme the color scheme used for the points. Possible color
##' schemes are: \code{\link{rainbow}}, \code{\link{heat.colors}},
##' \code{\link{topo.colors}}, \code{\link{terrain.colors}}.
##' @param point.alpha a \code{numerical value} indicating the transparency of
##' the point colors on a scale from 0 (transparent) to 1 (opaque).
##' @param cex \code{numerical value} for the size of the points.
##' @param line.col An character expression (any of \code{\link{colors}} or
##' hexadecimal notation), or numeric indicating the color of the line
##' connecting the point locations.
##' @return This function returns no data. It creates a .kml file in the in the
##' defined path.
##' @author Simeon Lisovski and Michael U. Kemp
##' @examples
##' \donttest{
##' data(hoopoe2)
##'   hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
##'   hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
##' filter <- distanceFilter(hoopoe2,distance=30)
##' ## takes time
##' ## trip2kml("trip.kml", hoopoe2$tFirst[filter], hoopoe2$tSecond[filter], hoopoe2$type[filter],
##' ##		degElevation=-6, col.scheme="heat.colors", cex=0.7,
##' ##		line.col="goldenrod")
##'}
##' @importFrom grDevices rgb
##' @export trip2kml
trip2kml <- function(file, tFirst, tSecond, type, degElevation, col.scheme="heat.colors", point.alpha=0.7, cex=1, line.col="goldenrod")
{
  if((length(tFirst)+length(type))!=(length(tSecond)+length(type))) stop("tFirst, tSecond and type must have the same length.")
  
  coord   <- coord(tFirst,tSecond,type,degElevation,note=F)
  index <- !is.na(coord[,2])
  datetime <- as.POSIXct(strptime(paste(ifelse(type==1,substring(tFirst,1,10),substring(tSecond,1,10)),
                                        " ",ifelse(type==1,"12:00:00","00:00:00"),sep=""),format="%Y-%m-%d %H:%M:%S"),"UTC")
  
  coord   <- coord[index,]
  longitude <- coord[,1]
  latitude <- coord[,2]
  
  date <- unlist(strsplit(as.character(datetime[index]), split = " "))[seq(1,
                                                                           ((length(datetime[index]) * 2) - 1), by = 2)]
  time <- unlist(strsplit(as.character(datetime[index]), split = " "))[seq(2,
                                                                           ((length(datetime[index]) * 2)), by = 2)]
  
  if(length(!is.na(coord[,2]))<1) stop("Calculation of coordinates results in zero spatial information.")
  
  if ((col.scheme%in% c("rainbow", "heat.colors", "terrain.colors", "topo.colors",
                        "cm.colors"))==F) stop("The col.scheme has been misspecified.")
  
  seq   <- seq(as.POSIXct(datetime[1]),as.POSIXct(datetime[length(datetime)]),by=12*60*60)
  index2<- ifelse(!is.na(merge(data.frame(d=datetime[index],t=TRUE),data.frame(d=seq,t=FALSE),by="d",all.y=T)[,2]),TRUE,FALSE)
  
  usable.colors <- strsplit(eval(parse(text = paste(col.scheme,
                                                    "(ifelse(length(index2) < 1000, length(index2), 1000), alpha=point.alpha)",
                                                    sep = ""))), split = "")[index2]
  
  usable.line.color <- strsplit(rgb(col2rgb(line.col)[1,
                                                      1], col2rgb(line.col)[2, 1], col2rgb(line.col)[3,
                                                                                                     1], col2rgb(line.col, alpha = 1)[4, 1], maxColorValue = 255),
                                split = "")
  
  date <- unlist(strsplit(as.character(datetime), split = " "))[seq(1,
                                                                    ((length(datetime) * 2) - 1), by = 2)]
  time <- unlist(strsplit(as.character(datetime), split = " "))[seq(2,
                                                                    ((length(datetime) * 2)), by = 2)]
  scaling.parameter <- rep(cex, length(latitude))
  
  data.variables <- NULL
  filename <- file
  write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", filename)
  write("<kml xmlns=\"http://www.opengis.net/kml/2.2\">", filename,
        append = TRUE)
  write("<Document>", filename, append = TRUE)
  write(paste("<name>", filename, "</name>", sep = " "),
        filename, append = TRUE)
  
  write("  <open>1</open>", filename, append = TRUE)
  write("\t<description>", filename, append = TRUE)
  write("\t  <![CDATA[Generated using <a href=\"http://simeonlisovski.wordpress.com/geolight\">GeoLight</a>]]>",
        filename, append = TRUE)
  write("\t</description>", filename, append = TRUE)
  write("<Folder>", filename, append = TRUE)
  write("  <name>Points</name>", filename, append = TRUE)
  write("<open>0</open>", filename, append = TRUE)
  for (i in 1:length(latitude)) {
    write("<Placemark id='point'>", filename, append = TRUE)
    write(paste("<name>", as.character(as.Date(datetime[i])), "</name>", sep = ""),
          filename, append = TRUE)
    write("  <TimeSpan>", filename, append = TRUE)
    write(paste("    <begin>", date[i], "T", time[i], "Z</begin>",
                sep = ""), filename, append = TRUE)
    write(paste("    <end>", date[ifelse(i == length(latitude),
                                         i, i + 1)], "T", time[ifelse(i == length(latitude),
                                                                      i, i + 1)], "Z</end>", sep = ""), filename, append = TRUE)
    write("  </TimeSpan>", filename, append = TRUE)
    write("<visibility>1</visibility>", filename, append = TRUE)
    write("<description>", filename, append = TRUE)
    write(paste("<![CDATA[<TABLE border='1'><TR><TD><B>Variable</B></TD><TD><B>Value</B></TD></TR><TR><TD>Date/Time</TD><TD>",
                datetime[i], "</TD></TR><TR><TD>lat long</TD><TD>",
                paste(latitude[i], longitude[i], sep = " "),
                "</TABLE>]]>", sep = "", collapse = ""), filename,
          append = TRUE)
    write("</description>", filename, append = TRUE)
    write("\t<Style>", filename, append = TRUE)
    write("\t<IconStyle>", filename, append = TRUE)
    write(paste("\t\t<color>", paste(noquote(usable.colors[[i]][c(8,
                                                                  9, 6, 7, 4, 5, 2, 3)]), collapse = ""), "</color>",
                sep = ""), filename, append = TRUE)
    write(paste("  <scale>", scaling.parameter[i], "</scale>",
                sep = ""), filename, append = TRUE)
    write("\t<Icon>", filename, append = TRUE)
    write("\t\t<href>http://maps.google.com/mapfiles/kml/pal2/icon26.png</href>",
          filename, append = TRUE)
    write("\t</Icon>", filename, append = TRUE)
    write("\t</IconStyle>", filename, append = TRUE)
    write("\t</Style>", filename, append = TRUE)
    write("\t<Point>", filename, append = TRUE)
    write(paste("\t<altitudeMode>", "relativeToGround", "</altitudeMode>",
                sep = ""), filename, append = TRUE)
    write("<tesselate>1</tesselate>", filename, append = TRUE)
    write("<extrude>1</extrude>", filename, append = TRUE)
    write(paste("\t  <coordinates>", longitude[i], ",", latitude[i],
                ",", 1,
                "</coordinates>", sep = ""), filename, append = TRUE)
    write("\t</Point>", filename, append = TRUE)
    write(" </Placemark>", filename, append = TRUE)
  }
  
  write("</Folder>", filename, append = TRUE)
  write("<Placemark>", filename, append = TRUE)
  write("  <name>Line Path</name>", filename, append = TRUE)
  write("  <Style>", filename, append = TRUE)
  write("    <LineStyle>", filename, append = TRUE)
  write(paste("\t<color>", paste(noquote(usable.line.color[[1]][c(8,
                                                                  9, 6, 7, 4, 5, 2, 3)]), collapse = ""), "</color>", sep = ""),
        filename, append = TRUE)
  
  
  write(paste("      <width>1</width>", sep = ""), filename,
        append = TRUE)
  write("    </LineStyle>", filename, append = TRUE)
  write("  </Style>", filename, append = TRUE)
  write("  <LineString>", filename, append = TRUE)
  write("    <extrude>0</extrude>", filename, append = TRUE)
  write("    <tessellate>1</tessellate>", filename, append = TRUE)
  write(paste("\t<altitudeMode>clampToGround</altitudeMode>",
              sep = ""), filename, append = TRUE)
  write(paste("     <coordinates>", noquote(paste(longitude,
                                                  ",", latitude, sep = "", collapse = " ")), "</coordinates>",
              sep = ""), filename, append = TRUE)
  write("    </LineString>", filename, append = TRUE)
  write("</Placemark>", filename, append = TRUE)
  write("</Document>", filename, append = TRUE)
  write("</kml>", filename, append = TRUE)
}

#' Draw the positions and the trip on a map
#'
#' Draw a map (from the \code{R} Package \code{maps}) with calculated positions
#' connected by a line
#'
#'
#' @param crds a \code{SpatialPoints} or \code{matrix} object, containing x
#' and y coordinates (in that order).
#' @param equinox logical; if \code{TRUE}, the equinox period(s) is shown as a
#' broken blue line.
#' @param ... Arguments to be passed to methods, such as graphical parameters (see par).
#' @param legend \code{logical}; if \code{TRUE}, a legend will be added to the plot.
#' @param xlim Longitude limits of map.
#' @param ylim Lattitude limits of map.
#' @author Simeon Lisovski
#' @examples
#'
#' data(hoopoe2)
#' hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
#' hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
#' crds <- coord(hoopoe2, degElevation = -6)
#' tripMap(crds, xlim = c(-20,20), ylim = c(0,60), main="hoopoe2")
#'
#' @importFrom ggplot2 ggtitle xlab ylab
#' @importFrom sf st_make_valid st_as_sf st_set_crs st_sf st_sfc st_linestring st_coordinates st_point st_combine
#' @importFrom dplyr summarize  %>% 
#' @export tripMap
tripMap <- function(crds, equinox=TRUE, xlim = NULL, ylim = NULL, legend = TRUE, ...) {
  
  args <- list(...)
  
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  dat <- crds %>% 
    as_tibble() %>% 
    filter(!is.nan(.data$lat))
  
  if(is.null(xlim) | is.null(ylim)) {
    world_sub <- world %>% 
      st_make_valid() %>% 
      st_intersection(st_bbox(dat %>% 
                                st_as_sf(coords = c('lon', 'lat'), crs = 4326)) %>%
                        st_as_sfc()) %>% 
      suppressWarnings() %>% suppressMessages()
  } else {
    world_sub <- world %>% st_make_valid() %>%
      st_intersection(st_bbox(c(xmin = xlim[1], 
                                xmax = xlim[2],
                                ymin = ylim[1],
                                ymax = ylim[2]), 
                              crs = 4326) %>%
                        st_as_sfc()) %>% 
      suppressWarnings() %>% 
      suppressMessages()
  }
  
  # Base map
  base <- ggplot() +
    theme_bw() +
    labs(x = ifelse(sum(names(args) %in% "xlab") == 1, args$xlab, "Longitude"),
         y = ifelse(sum(names(args) %in% "ylab") == 1, args$ylab, "Latitude"),
         title = ifelse(sum(names(args) %in% "main") == 1, args$main, "")) +
    geom_sf(data = world_sub, fill = "lightgray", color = "white") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) 
  
  
  # crds to sf
  points_sf <- dat %>%
    na.omit() %>% 
    st_as_sf(coords = c("lon","lat")) %>% 
    st_set_crs(4326) 
  
  # crosses
  pointPlot <- base +
    geom_sf(data = points_sf, mapping = aes(geometry = .data$geometry),
            size = ifelse(any(names(args) == "size"),
                          args$size, 1),
            shape = ifelse(any(names(args) == "shape"), args$shape, 3),
            show.legend = ifelse(any(names(args) == "show.legend"), args$show.legend, FALSE)) 
  
  # lines between crosses
  line_sf <- points_sf %>% 
    summarize() %>% 
    st_cast("LINESTRING") %>% 
    suppressMessages()
  
  linePlot <- pointPlot +
    geom_sf(data = line_sf, mapping = aes(geometry = .data$geometry),
            color = ifelse(sum(names(args)%in%"col")==1, args$col, "grey10"),
            linewidth = ifelse(sum(names(args)%in%"linewidth")==1, args$linewidth, 0.5),
            linetype = ifelse(sum(names(args)%in%"linetype")==1, args$linetype, 1)) 
  
  if(equinox){
    
    dat_eq <- apply(cbind(which(is.nan(crds[,2]) & !is.na(c(NaN, crds[-nrow(crds),2])))-1,
                          which(is.nan(crds[,2]) & !is.na(c(crds[-1,2], NaN)))+1), 1, 
                    function(x) st_as_sf(crds[c(x[1], x[2]),] %>% 
                                           as_tibble(), coords = c('lon', 'lat'), crs = 4326) %>%
                      st_combine() %>% st_cast("LINESTRING")) %>% Reduce("c",.)
    
    # Plot the ggplot with sf geometries
    equiPlot <- linePlot +
      geom_sf(data = dat_eq, mapping = aes(geometry = .data$geometry),
              color = "blue",
              linewidth = ifelse(sum(names(args)%in%"linewidth")==1, args$linewidth, 0.5), 
              linetype = ifelse(sum(names(args)%in%"linetype")==1, args$linetype, 1)) +
      labs(color = "Legend") 
    
  }
  # Check if legend should be added
  if(legend) {
    
    max_x <- max(st_coordinates(points_sf)[, "X"]) - 3
    min_y <- min(st_coordinates(points_sf)[, "Y"])
    
    independent_point <- st_sf(geometry = st_sfc(st_point(c(max_x - 10, min_y + 5)))) %>% 
      st_set_crs(4326)
    
    if(equinox) {
      
      legendPlot <- equiPlot + 
        geom_text(data = independent_point, 
                  aes(geometry = .data$geometry, label = "+ Position"), 
                  stat = "sf_coordinates", 
                  vjust = 0, 
                  hjust = 0) +
        geom_text(data = independent_point, 
                  aes(geometry = .data$geometry, label = "_____ Trip"), #used to be ——— but is not ASCII character 
                  stat = "sf_coordinates", 
                  vjust = 1.2, 
                  hjust = 0) +
        geom_text(data = independent_point, 
                  aes(geometry = .data$geometry, label = "_____ Equinox"), #used to be ——— but is not ASCII character
                  col = "blue", 
                  stat = "sf_coordinates", 
                  vjust = 2.4, 
                  hjust = 0) 
      
    } else {
      
      legendPlot <- linePlot + 
        geom_text(data = independent_point, 
                  aes(geometry = .data$geometry, label = "+ Position"), 
                  stat = "sf_coordinates", 
                  vjust = 0, 
                  hjust = 0) +
        geom_text(data = independent_point, 
                  aes(geometry = .data$geometry, label = "____ Trip"), 
                  stat = "sf_coordinates", 
                  vjust = 1.2, 
                  hjust = 0) 
      
    }
    print(legendPlot) %>% 
      suppressMessages() %>% 
      suppressWarnings()
    
  } else { #if not plot without legend
    if(ggplot) {
      if(equinox) {
        print(equiPlot) %>% 
          suppressWarnings() %>% 
          suppressMessages()
      } else {
        print(linePlot) %>% 
          suppressWarnings() %>% 
          suppressMessages()
      }
    }
  }
}

#' Calculate twilight events (sunrise/sunset) from light intensity measurements
#' over time
#'
#' Defines twilight events (sunrise/sunset) at times when the light intensity
#' measurements (\emph{light}) pass the defined light intensity threshold. An
#' interactive plot can be drawn to assess the calculations and improve e.g.
#' select only the realistic events.
#'
#'
#' @param datetime date and time of light intensity measurements e.g.
#' 2008-12-01 08:30 "UTC" (see:
#' \code{\link{as.POSIXct}},\link[=Sys.timezone]{time zones}).
#' @param light \code{numerical} value of the light intensity (usually
#' arbitrary units).
#' @param preSelection codelogical, if TRUE a pre selection of all calculated
#' twiligth events will be offered within the interactive process (ask=TRUE).
#' @param LightThreshold the light intensity threshold for the twilight event
#' calibration. If \code{Default}, it will be set slightly above (3 units) the
#' baseline level (measurement during the night).
#' @param maxLight if the geolocator record the maximum light value of a
#' certain time span, give the interval of maximum recordings in minutes (e.g.
#' 5).
#' @param ask \code{logical}, if TRUE the interactive plot will start after the
#' calculation.
#' @param nsee number of points to plot per screen.
#' @param allTwilights \code{logical}, if TRUE the function returns a list with
#' two tables
#' @return if allTwilights=FALSE, a \code{data frame}. Each row contains two
#' subsequent twilight events (\emph{tFirst, tSecond}) and \emph{type} defining
#' wether \emph{tFirst} refers to sunrise (1) or sunset (2). If
#' allTwilights=TRUE, a \code{list} with the data frame described in the
#' previous sentence and a data frame with all light intensities and a column
#' describing whether each row refers to sunrise (1), sunset (2) or to none of
#' these categories (0).
#' @note Depending on shading during light intensity measurements (e.g. due to
#' vegetation, weather, etc., see Lisovski et \emph{al.} 2012) the light
#' intensities may pass the light intensity threshold several times during the
#' day, resulting false sunrises and sunsets. It is highly recommended to check
#' the derived events visually (\code{ask=TRUE}).Twilight events can be deleted
#' and undeleted by clicking the (first) mouse button at the particular
#' position in the graph. The second mouse buttom (or esc) moves the time
#' series forward. Note, that a backward option is not included.
#' @author Simeon Lisovski
#' @export twilightCalc
#' @importFrom grDevices graphics.off
#' @importFrom graphics abline axis identify legend plot points 
twilightCalc <- function(datetime, light, LightThreshold=TRUE, preSelection=TRUE, maxLight=NULL, ask=TRUE, nsee=500, allTwilights=FALSE)
{
  if(class(datetime)[1]!="POSIXct") {
    stop(sprintf("datetime need to be provided as POSIXct class object."), call. = F)
    
  } else {
    bas <- data.frame(datetime=as.POSIXct(as.character(datetime),"UTC"),light)
  }
  
  
  if (is.numeric(LightThreshold))
  {
    LightThreshold <- as.numeric(LightThreshold)
    min <- min(bas$light)
  } else {
    # Basic level
    r <- as.data.frame(table(bas$light))
    nr <- as.numeric(which.max(r$Freq[as.numeric(r[,1])<mean(bas$light)]))
    LightThreshold <- (as.numeric(as.character(r[nr,1])))+3
  }
  
  out <- i.preSelection(bas$datetime,bas$light, LightThreshold)[,-1]
  
  if(!preSelection) out$mod <- 0
  
  if(ask)
  {
    n   <- nrow(bas)
    nn  <- n%/%nsee + 1
    cutsub <- cut(1:n, nn)
    picks <- NULL
    
    for(i in 1:nn){
      sub <- cutsub == levels(cutsub)[i]
      
      repeat{
        plot(bas[sub,1],bas[sub,2],type="o",cex=0.6,pch=20,ylab="Light intensity",xaxs="i",xaxt="n",xlab="",
             main=paste(as.Date(min(bas[sub,1]))," to ", as.Date(max(bas[sub,1]))," (end: ",as.Date(max(bas[,1])),")",sep=""))
        abline(h=LightThreshold,col="blue",lty=2)
        abline(v=out[out$mod==0,1],col="orange",lty=ifelse(out[out$mod==0,2]==1,1,2))
        points(out[,1],rep(LightThreshold,nrow(out[,])),col=ifelse(out$mod==0,"orange","grey"),pch=20,cex=0.8)
        
        axis(1,at=out[seq(from=1,to=nrow(out),length.out=(nrow(out)%/%2)),1],
             labels=substring(as.character(out[seq(from=1,to=nrow(out),length.out=(nrow(out)%/%2)),1]),6,16),cex=0.7)
        
        legend("topright",lty=c(3,1,2),lwd=c(1.3,2,2),col=c("blue",rep("orange",2)),c("Light\nThreshold","sunrise","sunset"),cex=1,bg="white")
        
        nr <- identify(out[,1],rep(LightThreshold,nrow(out)),n=1,plot=F)
        ifelse(length(nr)>0,ifelse(out$mod[nr]==0,out$mod[nr]<-1,out$mod[nr]<-0),break)
      }
    }
    
    cat("Thank you!\n\n")
    graphics.off()
  }
  
  results <- list()
  
  
  out <- subset(out,out$mod==0)[,-3]
  
  raw <- data.frame(datetime=c(as.POSIXct(datetime,"UTC"),as.POSIXct(out$datetime,"UTC")),
                    light=c(light,rep(LightThreshold,nrow(out))),type=c(rep(0,length(datetime)),out$type))
  
  raw <- raw[order(raw$type),]
  raw <- raw[-which(duplicated(as.character(raw$datetime),fromLast=T)),]
  raw <- raw[order(raw$datetime),]
  
  results$allTwilights <- raw
  
  opt <- data.frame(tFirst=as.POSIXct("1900-01-01 01:01","UTC"),tSecond=as.POSIXct("1900-01-01 01:01","UTC"),type=0)
  row <- 1
  for (k in 1:(nrow(out)-1))
  {
    if (as.numeric(difftime(out[k,1],out[k+1,1]))< 24 & out[k,1] != out[k+1,1])
    {
      opt[row,1] <- out[k,1]
      opt[row,2] <- out[k+1,1]
      opt[row,3] <- out$type[k]
      
      row <- row+1
    }
  }
  
  if(is.numeric(maxLight))
  {
    opt$tFirst[opt$type==2] <- opt$tFirst[opt$type==2] - (maxLight*60)
    opt$tSecond[opt$type==1] <- opt$tSecond[opt$type==1] - (maxLight*60)
  }
  
  if(allTwilights) {
    results$consecTwilights <- opt
    return(results)
  } else {
    return (opt)}
}

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
#'  calib2$tFirst <- as.POSIXct(calib2$tFirst, tz = "GMT")
#'  calib2$tSecond <- as.POSIXct(calib2$tSecond, tz = "GMT")
#' getElevation(calib2, known.coord = c(8,47.01))
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
#'  hoopoe2$tFirst <- as.POSIXct(hoopoe2$tFirst, tz = "GMT")
#'  hoopoe2$tSecond <- as.POSIXct(hoopoe2$tSecond, tz = "GMT")
#' coord <- coord(hoopoe2, degElevation=-6)
#' ## plot in a map using package maps
#' # par(oma=c(5,0,0,0))
#' # map(xlim=c(-20,40),ylim=c(-10,60),interior=F,col="darkgrey")
#' # map(xlim=c(-20,40),ylim=c(-10,60),boundary=F,lty=2,col="darkgrey",add=T)
#' # mtext(c("Longitude (degrees)","Latitude (degrees)"),side=c(1,2),line=c(2.2,2.5),font=3)
#' # map.axes()
#' # points(coord,col="brown",cex=0.5,pch=20)
#'
NULL