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

