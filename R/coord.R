#' Threshold based geographical positioning
#' 
#' Calculate the latitude and longitude from two subsequent twilight events
#' (sunrise/sunset)
#' 
#' The format (date and time) of \emph{tFirst} and \emph{tSecond} has to be
#' "yyyy-mm-dd hh:mm" corresponding to Universal Time Zone UTC (see:
#' \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones}).
#' 
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
#' respectively
#' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
#' twilight"). Either a single value, a \code{vector} with the same length as
#' \code{tFirst}. If site=TRUE a \code{vector} with a sun elevation angle for
#' each site in numerical order of the sites (incl. 0 for movements).
#' @param site \code{logical},if TRUE the a sun elevation angle for each site
#' (incl. 0 dor movements) will be used for positioning
#' @param sites a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0.
#' @param note \code{logical}, if TRUE a notation of how many positions could
#' be calculated in proportion to the number of failures will be printed at the
#' end.
#' @return A matrix of coordinates in decimal degrees. First column are
#' longitudes, expressed in degrees east of Greenwich. Second column contains
#' the latitudes in degrees north the equator. If latitude can not be
#' calculated (e.g. during equinox, at high latitudes) the function will return
#' NA for the date.
#' @author Simeon Lisovski
#' @references Montenbruck, O. & Pfleger, T. (2000) Astronomy on the Personal
#' Computer. \emph{Springer}, Berlin.
#' @examples
#' 
#' data(hoopoe2)
#' attach(hoopoe2)
#' coord <- coord(tFirst,tSecond,type,degElevation=-6)
#' ## tripMap(coord,xlim=c(-20,20),ylim=c(5,50), main="hoopoe2")
#' 
#' @export coord
coord <- function(tFirst,tSecond,type,degElevation=-6,site=FALSE,sites,note=TRUE) {

		
	# if noon, RadHourAngle = 0, if midnight RadHourAngle = pi
	# --------------------------------------------------------
	RadHourAngle <- numeric(length(type))
	index1 <- type==1
	if (sum(index1)>0) RadHourAngle[index1] <- 0
	RadHourAngle[!index1] <- pi
	# --------------------------------------------------------
		
	# --------------------------------------------------------
	Break <- FALSE
	if(site==TRUE){
		if(levels(as.factor(sites))!=length(degElevation)){
			cat("Error: Length of degElevation and number of sites (incl. 0 for movement) are not equal.")
		Break <- TRUE
		break
		}
	degEl <- rep(degElevation[1],length(tFirst))
	for(i in 1:length(sites[site!=0])){
		degEl[sites==i] <- 	degElevation[i]
	}
	degElevation <- degEl
	}
	# --------------------------------------------------------
	
	if(Break==FALSE){	
	
		tFirst  <- as.POSIXct(tFirst, "UTC")
		tSecond <- as.POSIXct(tSecond,"UTC")
		
		tSunTransit <- as.character(tFirst + as.numeric(difftime(tSecond,tFirst,units="secs")/2))
					   
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

out <- matrix(c(degLongitude,degLatitude),length(degLongitude),2)
if(note==TRUE) cat(paste("Note: Out of ",nrow(out)," twilight pairs, the calculation of ", nrow(out) - nrow(na.omit(out))," latitudes failed ","(",floor(((nrow(out)-nrow(na.omit(out)))*100)/nrow(out))," %)",sep=""))
return(out)

}
}

