
#' The GeoLight Package
#'
#' This is a summary of all features of \bold{\code{GeoLight}}, a
#' \code{R}-package for analyzing light based geolocator data.
#'
#' @name GeoLight-package
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




#' Residency analysis using a changepoint model
#'
#' Function to discriminate between periods of residency and movement based on
#' consecutive sunrise and sunset data. The calculation is based on a
#' changepoint model (\bold{\pkg{R}} Package \code{\link{changepoint}}:
#' \code{\link{cpt.mean}}) to find multiple changepoints within the
#' data.
#'
#' The \code{cpt.mean} from the codeR Package \code{changepoint} is a
#' function to find a multiple changes in mean for data where no assumption is
#' made on their distribution. The value returned is the result of finding the
#' optimal location of up to Q changepoints (in this case as many as possible)
#' using the cumulative sums test statistic.
#'
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
#' respectively
#' @param quantile probability threshold for stationary site selection. Higher
#' values (above the defined quantile of all probabilities) will be considered
#' as changes in the behavior. If specified, \code{rise.prob} and/or
#' \code{set.prob} will not be considered.
#' @param rise.prob the probability threshold for \bold{sunrise}: greater or
#' equal values indicates changes in the behaviour of the individual. If,
#' \code{NA} sunrise will not be considered.
#' @param set.prob the probability threshold for \bold{sunset}: higher and
#' equal values indicates changes in the behaviour of the individual. If,
#' \code{NA} sunrise will not be considered.
#' @param days a threshold for the length of stationary period. Periods smaller
#' than "days" will not be considered as a residency period
#' @param plot logical, if \code{TRUE} a plot will be produced
#' @param summary logical, if \code{TRUE} a summary of the results will be
#' printed
#' @return A \code{list} with probabilities for \emph{sunrise} and
#' \emph{sunset} the user settings of the probabilities and the resulting
#' stationary periods given as a \code{vector}, with the residency sites as
#' positiv numbers in ascending order (0 indicate movement/migration).
#' @note The sunrise and/or sunset times shown in the graph (if
#' \code{plot=TRUE}) represent hours of the day. However if one or both of the
#' twilight events cross midnight during the recording period the values will
#' be transformed to avoid discontinuity.
#' @author Simeon Lisovski & Tamara Emmenegger
#' @seealso \code{\link{changepoint}}, \code{\link{cpt.mean}}
#' @references Taylor, Wayne A. (2000) Change-Point Analysis: A Powerful New
#' Tool For Detecting Changes.
#'
#' M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis.
#' \emph{Wiley}.
#'
#' Chen, J. and Gupta, A. K. (2000) Parametric statistical change point
#' analysis. \emph{Birkhauser}.
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' residency <- changeLight(tFirst,tSecond,type,quantile=0.9)
#'
#' @export changeLight
#' @importFrom changepoint cpt.mean
changeLight <- function(tFirst,tSecond,type,quantile=0.6,rise.prob=NA,set.prob=NA,days=5,plot=TRUE,summary=TRUE)
{
	# start: Sunrise and Sunset
	rtype <- rep(2,length(type)); rtype[type==2] <- 1
	twE1  <- data.frame(ev=c(as.character(tFirst),as.character(tSecond)),t=c(type,rtype))
	twE2  <- twE1[!duplicated(twE1$ev),]

	twE2$ev <- as.POSIXct(as.character(twE2$ev),"UTC")

	twE <- twE2[order(twE2$ev),]

	rise <- as.numeric(substring(twE$ev[twE$t==1],12,13)) +  as.numeric(substring(twE$ev[twE$t==1],15,16))/60
	set  <- as.numeric(substring(twE$ev[twE$t==2],12,13)) +  as.numeric(substring(twE$ev[twE$t==2],15,16))/60
	# end: Sunrise and Sunset

	cor.rise <- rep(NA, 24)
	for(i in 0:23){
		cor.rise[i+1] <- max(abs((c(rise[1],rise)+i)%%24 -
		            (c(rise,rise[length(rise)])+i)%%24),na.rm=T)
	}
	rise <- (rise + (which.min(round(cor.rise,2)))-1)%%24

	cor.set <- rep(NA, 24)
	for(i in 0:23){
		cor.set[i+1] <- max(abs((c(set[1], set)+i)%%24 -
		            (c(set, set[length(set)])+i)%%24),na.rm=T)
	}
	set <- (set + (which.min(round(cor.set,2)))-1)%%24


	# start: Change Point Model
	# max. possible Change Points (length(sunrise)/2)
	CPs1 <- cpt.mean(rise, Q=length(rise)/2, pen.value=0.001)
	CPs2 <- cpt.mean(set, Q=length(set)/2, pen.value=0.001)

N1 <- seq(1,length(rise))
N2 <- seq(1,length(set))

	tab1 <- merge(data.frame(N=N1,prob=NA),data.frame(N=CPs1$cps[1,],prob=CPs1$cps[2,]),by.x="N",by.y="N",all.x=T)[,-2]
		tab1[is.na(tab1[,2]),2] <- 0
	tab2 <- merge(data.frame(N=N2,prob=NA),data.frame(N=CPs2$cps[1,],prob=CPs2$cps[2,]),by.x="N",by.y="N",all.x=T)[,-2]
		tab2[is.na(tab2[,2]),2] <- 0
	# end: Change Point Model

# quantile calculation
if(is.numeric(quantile))
{
rise.prob<-as.numeric(round(as.numeric(quantile(tab1[tab1[,2]!=0,2],probs=quantile,na.rm=TRUE)),digits=5))
set.prob<-as.numeric(round(as.numeric(quantile(tab2[tab2[,2]!=0,2],probs=quantile,na.rm=TRUE)),digits=5))
}

if(plot){
	layout(matrix(c(4,1,2,3),nrow=4,byrow=T),heights=c(0.5,1,0.5,0.5))

	par(mar=c(2,4.5,2,5),cex.lab=1.5,cex.axis=1.5,bty="o")
	plot(twE$ev[twE$t==1],rise,type="o",cex=0.2,col="red",ylab="Sunrise (red)",xlim=c(min(twE[,1]),max(twE[,1])),xaxt="n")
	par(new=T)
	plot(twE$ev[twE$t==2],set,type="o",cex=0.2,col="blue",xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(min(twE[,1]),max(twE[,1])))
	axis(4)
	mtext("Sunset (blue)",4,line=2.7,cex=1)
	axis(1,at=seq(min(twE[,1]),max(twE[,1]),by=(10*24*60*60)),labels=F)
	axis(1,at=seq(min(twE[,1]),max(twE[,1]),by=(30*24*60*60)),lwd.ticks=2,labels=as.Date(seq(min(twE[,1]),max(twE[,1]),by=(30*24*60*60))),cex.axis=1)


	par(mar=c(1.5,4.5,0.8,5),bty="n")
	plot(rep(twE$ev[twE$t==1][1],2),c(0,tab1[1,2]),type="l",lwd=3,col="red",ylab="",xaxt="n",xlim=c(min(twE[,1]),max(twE[,1])),ylim=c(0,max(na.omit(c(tab1[,2],tab2[,2])))))
	for(i in 2:nrow(tab1)){
		lines(rep(twE$ev[twE$t==1][i],2),c(0,tab1[i,2]),lwd=3,col="red")
	}

	if(is.numeric(rise.prob)){abline(h=rise.prob,lty=2)}
	par(mar=c(1.5,4.5,0.8,5),bty="n")

	plot(rep(twE$ev[twE$t==2][1],2),c(0,tab2[1,2]),type="l",lwd=3,col="blue",ylab="",xaxt="n",xlim=c(min(twE[,1]),max(twE[,1])),ylim=c(0,max(na.omit(c(tab1[,2],tab2[,2])))))
	for(i in 2:nrow(tab1)){
		lines(rep(twE$ev[twE$t==2][i],2),c(0,tab2[i,2]),lwd=3,col="blue")
	}
	if(is.numeric(set.prob)){abline(h=set.prob,lty=2)}
	mtext("Probability of change",side=2,at=max(na.omit(c(tab1[,2],tab2[,2]))),line=3)
}


	out <- list()
	out$riseProb <- tab1
		out$riseProb[out$riseProb==0] <- NA
	out$setProb  <- tab2
		out$setProb[out$setProb==0] <- NA

	# site calculation
	r <- out$riseProb[,2]
	s <- out$setProb[,2]

	if(is.na(rise.prob)) r <- rep(NA,nrow(out$riseProb))
	if(is.na(set.prob))  s <- rep(NA,nrow(out$setProb))



	comp <- data.frame(date=twE$ev, t= twE$t, prob=NA,set.probs=NA)
			 comp[twE$t==1,3] <- r ; comp[twE$t==2,3] <- s
			 comp[twE$t==1,4] <- rise.prob ; comp[twE$t==2,4] <- set.prob

	comp$filt <- FALSE
	comp$filt[comp[,3]>=comp$set.prob] <- TRUE


	# start: Site selection procedure
	site <- rep(0,nrow(comp))
	date <- comp$date
	if(!comp$filt[1]) site[1] <- 1
	for(i in 2:length(date)){

	if(!comp$filt[i] & site[i-1]!=0) site[i] <- as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])
	if(!comp$filt[i] & site[i-1]==0) site[i] <- as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])+1
	if(comp$filt[i]  & abs(as.numeric(
						difftime(
								min(date[site==as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])]),
								max(date[site==as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])]),
								units="days"
								)
								))
		<days){
		site[site==as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])] <- 0
		}
	}
	# end: Site selection procedure

	out$rise.prob <- rise.prob
	out$set.prob <- set.prob

# Schedule
	diff1 <- c(site,site[length(site)]) != c(0,site)
	a <- which(diff1)[c(T,F)]
	b <- which(diff1)[c(F,T)]

rows <- sort(c(
				(if(a[1]==1) c(a[1],a[-1]-1) else a-1)
				,b,length(c(out$riseProb[,2],out$setProb[,2]))))

date <- comp$date[rows]

	sc <- data.frame(site=c(letters[1:(length(date)/2)]),
		  start=as.POSIXct(date[c(T,F)],tz="UTC"),
		  end=as.POSIXct(date[c(F,T)],tz="UTC")
		  )

index <- sc$start[-1]==sc$end[-nrow(sc)]; sc$start[-1][index] <- (comp$date[rows+1][c(T,F)][-1])[index]
# end Schedule

# site translation to midnoon data
tFirst <- as.POSIXct(as.character(tFirst),"UTC")
tSecond<- as.POSIXct(as.character(tSecond),"UTC")

midnoon <- tFirst + (tSecond-tFirst)/2
mdSite  <- rep(0,length(midnoon))

for(i in 1:nrow(sc)){
	mdSite[midnoon>sc$start[i] & midnoon<sc$end[i]] <- i
}

	out$site <- mdSite


if(plot) {
	par(mar=c(1,4.5,1,5),bty="o")
		mig <- out$site
		mig[mig>0] <- 1
	plot(midnoon,mig,type="l",yaxt="n",ylab=NA,ylim=c(0,1.5))
		rect(sc$start,1.1,sc$end,1.4,lwd=0,col="grey")
	}

out$migTable <- data.frame(site=letters[1:nrow(sc)],
							arrival=as.Date(sc$start),
							departure=as.Date(sc$end),
							days=round(as.numeric(difftime(sc$end,sc$start,units="days")),0),
							P.start=round(comp$prob[rows][c(T,F)],4),
							P.end=round(comp$prob[rows][c(F,T)],4))


if(summary){i.sum.Cl(out)}

return(out)

}


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
	if(site){
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

	if(!Break){

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
if(note) cat(paste("Note: Out of ",nrow(out)," twilight pairs, the calculation of ", nrow(out) - nrow(na.omit(out))," latitudes failed ","(",floor(((nrow(out)-nrow(na.omit(out)))*100)/nrow(out))," %)",sep=""))
return(out)

}
}

#' Filter for unrealistic positions within a track
#'
#' The filter identifies unrealistic positions. The maximal distance per
#' hour/day can be set corresponding to the particular species.
#'
#' The filter currently depends to some extent on the calibration
#' (\code{degElevation}).
#'
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining tFirst as sunrise or sunset respectively
#' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
#' twilight")
#' @param distance the maximal distance in km per \code{units}. Distances above
#' will be considered as unrealistic.
#' @param units the time unite corresponding to the distance. Default is
#' "hour", alternative option is "day".
#' @return Logical \code{vector} matching positions that pass the filter.
#' @author Simeon Lisovski, Fraenzi Korner-Nievergelt
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' coord  <- coord(tFirst,tSecond,type)
#' filter <- distanceFilter(tFirst,tSecond,type,distance=30)
#' tripMap(coord[filter,],xlim=c(-20,20),ylim=c(0,60),main="hoopoe2 (filter)")
#'
#' @export distanceFilter
distanceFilter <- function(tFirst,tSecond,type,degElevation=-6,distance,units="hour") {

if(units=="days") units <- distance/24

tFirst <- as.POSIXct(as.character(tFirst),"UTC")
tSecond<- as.POSIXct(as.character(tSecond),"UTC")
tSunTransit <- tFirst + (tSecond-tFirst)/2
coord  <- coord(tFirst,tSecond,type,degElevation,note=FALSE)

coord[is.na(coord[,2]),2] <- 999

difft <- as.numeric(difftime(tSunTransit[-length(tSunTransit)],tSunTransit[-1],units="hours"))
diffs <- abs(as.numeric(i.loxodrom.dist(coord[-nrow(coord),1],coord[-nrow(coord),2],coord[-1,1],coord[-1,2]))/difft)

index <- rep(TRUE,length(tFirst))
index[diffs>distance] <- FALSE
index[coord[,2]==999] <- TRUE


cat(paste("Note: ",length(index[!index])," of ",length(index[coord[,2]!=999])," positions were filtered (",floor((length(index[!index])*100)/length(index[coord[,2]!=999]))," %)",sep=""))
return(index)
}


#' Calculate the appropriate sun elevation angle for known location
#'
#' Function to calculate the sun elevation angle for light measurements at a
#' known location
#'
#'
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining tFirst as sunrise or sunset respectively
#' @param known.coord a \code{SpatialPoint} or \code{matrix} object, containing
#' known x and y coordinates (in that order) for the selected measurement
#' period.
#' @param plot \code{logical}, if TRUE a plot will be produced.
#' @author Simeon Lisovski
#' @references Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt,
#' F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
#' precision affected by environmental factors. \emph{Methods in Ecology and
#' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
#' @examples
#'
#' data(calib2)
#' attach(calib2)
#' known.coord <- c(7.1,46.3)
#' getElevation(tFirst,tSecond,type,known.coord)
#'
#' @export getElevation
getElevation <- function(tFirst,tSecond,type,known.coord,plot=TRUE) {


 	table <- data.frame(tFirst=as.POSIXct(as.character(tFirst),"UTC"),tSecond = as.POSIXct(as.character(tSecond),"UTC"),type=as.numeric(type))

 coord <- known.coord
 tab <- table
 degElevation <- ifelse(coord[2]<0,-12,12)

 	lat0   <- coord(tab[,1],tab[,2],tab[,3],degElevation,note=F)[,2]
 	while(length(lat0[!is.na(lat0)])!=nrow(tab)){
 		degElevation <- degElevation - ifelse(coord[2]<0,-0.025,+0.025)
 		lat0  <- coord(tab[,1],tab[,2],tab[,3],degElevation,note=F)[,2]
 	}




 	lat1   <- coord(tab[,1],tab[,2],tab[,3],degElevation,note=F)[,2]
 	x0 <- i.loxodrom.dist(coord[1],coord[2],coord[1],median(lat1[!is.na(lat1)]))
  	degElevation <- degElevation - ifelse(coord[2]<0,-0.025,+0.025)
  	lat2   <- coord(tab[,1],tab[,2],tab[,3],degElevation,note=F)[,2]
  	x1 <- i.loxodrom.dist(coord[1],coord[2],coord[1],median(lat2[!is.na(lat2)]))

 	while(x0 > x1)
 		{
 			degElevation <- degElevation - ifelse(coord[2]<0,-0.025,+0.025)
 			if(degElevation==ifelse(coord[2]<0,10,-10)) break
 			x0 <- x1
 			latNew <- coord(tab[,1],tab[,2],tab[,3],degElevation,note=F)[,2]
 			x1 <- i.loxodrom.dist(coord[1],coord[2],coord[1],median(latNew[!is.na(latNew)]))
 		}

 if(plot)
 	{
 	Tw <- as.POSIXct(subset(tab,type=1,select=c(tFirst))$tFirst,"UTC")

 	SElev <- i.sunelevation(coord[1],coord[2],as.numeric(substring(Tw,1,4)),as.numeric(substring(Tw,6,7)),
         as.numeric(substring(Tw,9,10)),as.numeric(substring(Tw,12,13)),as.numeric(substring(Tw,15,16)),0)

 	layout(matrix(c(1,2,3,3),nrow=2,ncol=2,byrow=TRUE))
 	par(oma=c(0.2,0.2,2,0.2))
 	par(mar=c(6,4,6,1),bty="n",yaxt="n",xaxt="s")

 	plot(SElev[as.numeric(substring(Tw,12,13)) %in% 0:12],
         jitter(rep(1,length(SElev[as.numeric(substring(Tw,12,13)) %in% 0:12])),0.2),pch=20,
         cex=1,xlim=c(-10,max(SElev)+3),ylim=c(0.9,1.1),ylab="",xlab="",main="Sunrise",cex.main=1.1,font.main=3)
 	mtext("Light intensity threshold (jitter)",side=2,cex=1.1,font=6)
 	arrows(-9.8,1,-8.64,1,length=0.1)
 	abline(v=-6,lty=2,lwd=0.3)
 	abline(v=degElevation,lty=2,lwd=2,col="orange")


 	par(mar=c(6,1,6,3),bty="n",yaxt="n",xaxt="s")
 	plot(SElev[as.numeric(substring(Tw,12,13)) %in% 13:23],
       jitter(rep(1,length(SElev[as.numeric(substring(Tw,12,13)) %in% 13:23])),0.2),
       pch=1,cex=0.7,xlim=c(max(SElev)+3,-10),ylim=c(0.9,1.1),xlab="",main="Sunset",cex.main=1.1,font.main=3)
 	abline(v=-6,lty=2,lwd=0.3)
 	abline(v=degElevation,lty=2,lwd=2,col="orange")

 	legend("topleft",lty=c(2,2,2),lwd=c(0.3,2,2),col=c("black","transparent","orange"),
       c("- 6 degrees","",paste("getElevation\n",round(degElevation-0.025,3)," degrees",sep="")),bg="white",box.col="white",cex=.9)

 	mtext("Twilight times over sun elevation angles",line=0, adj=0.52, cex=1.5,col="black", outer=TRUE)

 	par(bty="o",mar=c(6,6,1,1),yaxt="s")
 	t <- seq(tab$tFirst[1],tab$tSecond[nrow(tab)],by=60)
 	plot(t,i.sunelevation(coord[1],coord[2],as.POSIXlt(t)$year+1900,as.POSIXlt(t)$mo+1,as.POSIXlt(t)$mday,as.POSIXlt(t)$hour,as.POSIXlt(t)$min,as.POSIXlt(t)$sec),type="l",
 		xaxt="n",xlab="",ylab="Sun elevation angle")
 	abline(h=degElevation,lty=2,col="grey80")
 	points(c(tab$tFirst,tab$tSecond),rep(degElevation,(nrow(tab)*2)),pch="+",col="darkgreen",cex=2)
	t2 <- seq(tab$tFirst[1],tab$tSecond[nrow(tab)],by=2*24*60*60)
	axis(1,at=t2,labels=as.character(as.Date(t2)))
 	}

 return(degElevation - 0.025)
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
#' @export glfTrans
glfTrans <-
function(file="/path/file.glf") {


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
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
#' respectively
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0
#' @param start.angle a single sun elevation angle. The combined process of
#' checking for minimal variance in resulting latitude, which is the initial
#' value for the sun elevation angle in the iterative process of identifying
#' the latitudes with the least variance
#' @param distanceFilter logical, if TRUE the \code{\link{distanceFilter}} will
#' be used to filter unrealistic positions
#' @param distance if \code{distanceFilter} is set \code{TRUE} a threshold
#' distance in km has to be set (see: \code{\link{distanceFilter}})
#' @param plot logical, if TRUE the function will give a plot with all relevant
#' information
#' @return A \code{vector} of sun elevation angles corresponding to the
#' Hill-Ekstrom calibration for each defined period.
#' @author Simeon Lisovski
#' @references Ekstrom, P.A. (2004) An advance in geolocation by light.
#' \emph{Memoirs of the National Institute of Polar Research}, Special Issue,
#' \bold{58}, 210-226.
#'
#' Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step:
#' Latitude. In: \emph{Electronic Tagging and Tracking in Marine Fisheries}
#' (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The
#' Netherlands.
#'
#' Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F.,
#' Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and
#' precision affected by environmental factors. \emph{Methods in Ecology and
#' Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' residency <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)
#' HillEkstromCalib(tFirst,tSecond,type,residency$site,-6)
#'
#' @export HillEkstromCalib
HillEkstromCalib <-
function(tFirst,tSecond,type,site,start.angle=-6,distanceFilter=FALSE,distance,plot=TRUE) {

	tFirst <- as.POSIXct(as.character(tFirst),"UTC")
	tSecond <- as.POSIXct(as.character(tSecond),"UTC")

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
	t0 <- var(na.omit(coord(tFirst[site==j],tSecond[site==j],type[site==j],start-((b*0.1)-0.1),note=F)[,2]))
	t1 <- var(na.omit(coord(tFirst[site==j],tSecond[site==j],type[site==j],start-(b*0.1),note=F)[,2]))
	t2 <- var(na.omit(coord(tFirst[site==j],tSecond[site==j],type[site==j],start-((b*0.1)+0.1),note=F)[,2]))
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
	f0 <- var(na.omit(coord(tFirst[site==j],tSecond[site==j],type[site==j],start+((b*0.1)-0.1),note=F)[,2]))
	f1 <- var(na.omit(coord(tFirst[site==j],tSecond[site==j],type[site==j],start+(b*0.1),note=F)[,2]))
	f2 <- var(na.omit(coord(tFirst[site==j],tSecond[site==j],type[site==j],start+((b*0.1)+0.1),note=F)[,2]))
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


layout(matrix(seq(1,sites*2),ncol=2,nrow=sites,byrow=T))
par(oma=c(3,5,6,2))

for(j in 1:sites){

	if(is.na(HECalib[j])){

	par(mar=c(2,2,2,2),bty="n")
	plot(0,0,cex=0,pch=20,col="white",ylab="",xlab="",xaxt="n",yaxt="n")
	text(0,0,"NA",cex=2)
	par(mar=c(2,2,2,2),bty="n")
	plot(0,0,cex=0,pch=20,col="white",ylab="",xlab="",xaxt="n",yaxt="n")

	} else {

	angles <- c(seq(HECalib[j]-2,HECalib[j]+2,0.2))

	latM <- matrix(ncol=length(angles),nrow=length(tFirst[site==j]))

	for(i in 1:ncol(latM)){
	latM[,i] <- coord(tFirst[site==j],tSecond[site==j],type[site==j],c(angles[i]),note=F)[,2]
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
	par(mar=c(2,2,2,2))
	plot(tFirst[site==j],latT[,1],ylim=c(min(na.omit(min)),max(na.omit(max))),type="o",cex=0.7,pch=20,col=colors[1],ylab="",xlab="")
	for(p in 2:ncol(latT)){
		   lines(tFirst[site==j],latM[,p],type="o",cex=0.7,pch=20,col=colors[p])
		   }
	lines(tFirst[site==j],coord(tFirst[site==j],tSecond[site==j],type[site==j],HECalib[j],note=F)[,2],col="tomato2",type="o",lwd=2,cex=1,pch=19)
	if(j==sites) mtext("Latitude",side=2,line=3)
	if(j==sites) mtext("Date",side=1,line=2.8)


	par(mar=c(2,2,2,2))
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


	return(HECalib)
}

i.JC2000 <-
function(jD) {

#--------------------------------------------------------------------------------------------------------
# jD: julian Date
#--------------------------------------------------------------------------------------------------------

options(digits=10)


	  	jC<- (jD - 2451545)/36525

return(jC)

			  }

i.deg <-
function(Rad) {

#------------------------------------------------------------
# Deg: 	The input angle in radian
#------------------------------------------------------------

options(digits=10)


     return(Rad * (180/pi))

}

i.frac <-
function(In) {

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

i.julianDate <-
function(year,month,day,hour,min) {

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

i.loxodrom.dist <-
function(x1, y1, x2, y2, epsilon=0.0001){
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


i.rad <-
function(Deg) {

#------------------------------------------------------------
# Deg: 	The input angle in degrees
#------------------------------------------------------------

options(digits=10)


   	return(Deg * (pi/180))

}

i.radDeclination <-
function(radEclipticalLength,radObliquity) {

#-------------------------------------------------------------------------------------------------------------------
# RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane.
# EadObliquity:
#-------------------------------------------------------------------------------------------------------------------

options(digits=10)

	dec <- asin(sin(radObliquity)*sin(radEclipticalLength))

return(dec)

}

i.radEclipticLongitude <-
function(jC) {

#-------------------------------------------------------------------------------------------------------------------
# jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00
#-------------------------------------------------------------------------------------------------------------------



	options(digits=10)

  	radMeanAnomaly <- 2*pi*i.frac(0.993133 + 99.997361*jC)
	EclipticLon    <- 2*pi*i.frac(0.7859452 + radMeanAnomaly/(2*pi) + (6893*sin(radMeanAnomaly) + 72*sin(2*radMeanAnomaly) + 6191.2*jC) / 1296000)

return(EclipticLon)

}

i.radGMST <-
function(jD,jD0,jC,jC0) {

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

i.radObliquity <-
function(jC) {

#--------------------------------------------------------------------------------------------------------
# jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
#--------------------------------------------------------------------------------------------------------

options(digits=10)


degObliquity <- 23.43929111 - (46.8150 + (0.00059 - 0.001813 *jC)*jC) *jC/3600
radObliquity <- i.rad(degObliquity)

return(radObliquity)

}

i.radRightAscension <-
function(RadEclipticalLength,RadObliquity) {

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

i.setToRange <-
function(Start,Stop,Angle) {

#-------------------------------------------------------------------------------------------------------------------
# Start:	 Minimal value of the range in degrees
# Stop: 	 Maximal value of the range in degrees
# Angle:	 The angle that should be fit to the range
#-------------------------------------------------------------------------------------------------------------------

options(digits=15)

		angle <- Angle
		range <- Stop - Start



			index1 <- angle >= Stop
			if (sum(index1)>0) angle[index1] <- angle[index1] - (floor((angle[index1]-Stop)/range)+1)*range

			index2 <- angle < Start
			if(sum(index2)>0) angle[index2] <- angle[index2]  + ceiling(abs(angle[index2] -Start)/range)*range

return(angle)

}

# Summery of the changeLight Function
# ---------------------------------------------------------------------------
# Author: Simeon Lisovski, Mai 2012
# ---------------------------------------------------------------------------

i.sum.Cl <- function(object) {

	if(sum(names(object)==c("riseProb","setProb","rise.prob","set.prob","site","migTable"))==6){
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

i.sunelevation <-
function(lon, lat, year, month, day, hour, min, sec){

#-------------------------------------------------------------------------------
# lon: longitude in decimal coordinates
# lat: latitude in decimal coordinates
# year: numeric, e.g. 2006 (GMT)
# month: numeric, e.g. 8  (GMT)
# day: numeric, e.g. 6   (GMT)
# hour:  numeric e.g. 6  (GMT)
# min: numeric, e.g. 0   (GMT)
# sec: numeric, e.g.     (GMT)
#-------------------------------------------------------------------------------

datetime<-paste(year,"-", month,"-", day, " ", hour, ":", min, ":", sec, sep="")
gmt<-as.POSIXct(strptime(datetime, "%Y-%m-%d %H:%M:%S"), "UTC")
n <- gmt - as.POSIXct(strptime("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S"), "UTC")

# mean ecliptical length of sun
L <- 280.46 + 0.9856474 * n
L <- as.numeric(L)

# Anomalie
g <- 357.528 + 0.9856003 * n
g <- as.numeric(g)

t.v <- floor(g/360)
g <- g - 360*t.v
g.rad <- g*pi/180

t.l <- floor(L/360)
L <- L - 360 * t.l
L.rad <- L*pi/180

# ecliptical length of sun
LAMBDA <- L + 1.915 * sin(g.rad) + 0.02*sin(2*g.rad)
LAMBDA.rad <- LAMBDA*pi/180

# coordinates of equator
epsilon <- 23.439 - 0.0000004 * n
epsilon.rad <- as.numeric(epsilon)*pi/180

alpha.rad <- atan(cos(epsilon.rad)*sin(LAMBDA.rad)/cos(LAMBDA.rad))


alpha.rad <- ifelse(cos(LAMBDA.rad)<0, alpha.rad+pi, alpha.rad)
alpha <- alpha.rad*180/pi

deklination.rad <- asin(sin(epsilon.rad) * sin(LAMBDA.rad))
deklination <- deklination.rad*180/pi

# angle h
tag<-paste(year,"-", month,"-", day, " 00:00:00", sep="")
JD0<-as.POSIXct(strptime(tag, "%Y-%m-%d %H:%M:%S"), "GMT")
JD0 <- JD0 - as.POSIXct(strptime("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S"), "GMT")
T0 <- JD0/36525

Time <- hour  + min/60  + sec/60/100
theta.Gh <- 6.697376 + 2400.05134 * T0 + 1.002738 * Time
theta.Gh <- as.numeric(theta.Gh)

t.d <- floor(theta.Gh/24)
theta.Gh <- theta.Gh-t.d*24

theta.G <- theta.Gh * 15

theta <- theta.G + lon      # Stundenwinkel des Fruhlingspunktes
tau <- theta-alpha    # Stundenwinkel
tau.rad <- tau/180*pi

# H?henwinkel h
h <- asin(cos(deklination.rad) * cos(tau.rad) * cos(lat/180*pi) + sin(deklination.rad) * sin(lat/180*pi))
h.grad <- h/pi*180

# correction because of refraction
R <- 1.02/(tan((h.grad+10.3/(h.grad+5.11))/180*pi))
hR.grad <- h.grad + R/60
return(hR.grad)
}

i.twilightEvents <- function(datetime, light, LightThreshold){

   df <- data.frame(datetime, light)

   ind1 <- which((df$light[-nrow(df)] < LightThreshold & df$light[-1] > LightThreshold) |
    			 (df$light[-nrow(df)] > LightThreshold & df$light[-1] < LightThreshold) |
  				  df$light[-nrow(df)] == LightThreshold)

   bas1 <- cbind(df[ind1,],df[ind1+1,])
  		  bas1 <- bas1[bas1[,2]!=bas1[,4],]

  x1 <- as.numeric(unclass(bas1[,1])); x2 <- as.numeric(unclass(bas1[,3]))
  y1 <- bas1[,2]; y2 <- bas1[,4]
  m <- (y2-y1)/(x2-x1)
  b <- y2-(m*x2)

  xnew <- (LightThreshold - b)/m
  type <- ifelse(bas1[,2]<bas1[,4],1,2)
  res  <- data.frame(datetime=as.POSIXct(xnew, origin="1970-01-01", tz="UTC"),type)

return(res)

}#' Filter to remove noise in light intensity measurements during the night
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

return(light)

}
#' Filter to remove outliers in defined sunrise and sunset times
#'
#' This filter defines outliers based on residuals from a local polynomial
#' regression fitting provcess (\code{\link{loess}}).
#'
#'
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
#' respectively
#' @param k a measure of how many interquartile ranges to take before saying
#' that a particular twilight event is an outlier
#' @param plot codelogical, if TRUE a plot indicating the filtered times will
#' be produced.
#' @return Logical \code{vector} matching positions that pass the filter.
#' @author Simeon Lisovski & Eldar Rakhimberdiev
#' @export loessFilter
loessFilter <- function(tFirst, tSecond, type, k=3, plot=TRUE){

	tw <- data.frame(datetime=as.POSIXct(c(tFirst,tSecond),"UTC"),type=c(type,ifelse(type==1,2,1)))
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
par(mfrow=c(2,1),mar=c(3,3,0.5,3),oma=c(2,2,0,0))
plot(dawn$datetime[dawn$type==1],dawn$hours[dawn$type==1],pch="+",cex=0.6,xlab="",ylab="",yaxt="n")
lines(dawn$datetime[!dawn$filter], predict(loess(dawn$hours[!dawn$filter]~as.numeric(dawn$datetime[!dawn$filter]),span=0.1)) , type="l")
points(dawn$datetime[dawn$filter],dawn$hours[dawn$filter],col="red",pch="+",cex=1)
axis(2,labels=F)
mtext("Sunrise",4,line=1.2)

plot(dusk$datetime[dusk$type==2],dusk$hours[dusk$type==2],pch="+",cex=0.6,xlab="",ylab="",yaxt="n")
lines(dusk$datetime[!dusk$filter], predict(loess(dusk$hours[!dusk$filter]~as.numeric(dusk$datetime[!dusk$filter]),span=0.1)), type="l")
points(dusk$datetime[dusk$filter],dusk$hours[dusk$filter],col="red",pch="+",cex=1)
axis(2,labels=F)
legend("bottomleft",c("Outside filter","Inside filter"),pch=c("+","+"),col=c("black","red"),
	   bty="n",cex=0.8)
mtext("Sunset",4,line=1.2)
mtext("Time",1,outer=T)
mtext("Sunrise/Sunset hours (rescaled)",2,outer=T)
}
all <- rbind(subset(dusk,filter),subset(dawn,filter))

filter <- rep(FALSE,length(tFirst))
	filter[tFirst%in%all$datetime | tSecond%in%all$datetime] <- TRUE

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
#' @export luxTrans
luxTrans <-
function(file="/path/file.lux") {

lux1 <- read.table(file,sep="\t",skip=21,col.names=c("datetime","time")) # read file
lux <- data.frame(datetime=as.POSIXct(strptime(lux1[,1],format="%d/%m/%Y %H:%M:%S",tz="UTC")),light=lux1[,2])

return(lux)
}

# Function to summerize the individual migration pattern (Package GeoLight)
# ---------------------------------------------------------------------------
# Author: Simeon Lisovski, November 2011
# ---------------------------------------------------------------------------



#' Summary of migration/movement pattern
#'
#' Function for making a data frame summarising residency and movement pattern.
#'
#'
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param site a \code{vector}, indicating the residency period of a particular
#' day (see output: \code{\link{changeLight}})
#' @return A \code{data.frame} with end and start date (yyyy-mm-dd hh:mm, UTC)
#' for each stationary period.
#' @author Simeon Lisovski
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' residency <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)
#' schedule(tFirst,tSecond,residency$site)
#'
#' @export schedule
schedule <- function(tFirst,tSecond,site){

	tFirst <- as.POSIXct(as.character(tFirst),"UTC")
	tSecond <- as.POSIXct(as.character(tSecond),"UTC")


	midnoon <- tFirst + (tSecond - tFirst)/2
	diff1 <- c(site,site[length(site)]) != c(0,site)
	a <- which(diff1)[c(T,F)]
	b <- which(diff1)[c(F,T)]

	rows <- sort(c(
					(if(a[1]==1) c(a[1],a[-1]-1) else a-1)
					,b,length(midnoon)))

	midnoon1 <- midnoon[rows]

	st <- data.frame(site=c(letters[1:(length(midnoon1)/2)]),
		  start=as.POSIXct(midnoon1[c(T,F)],tz="UTC"),
		  end=as.POSIXct(midnoon1[c(F,T)],tz="UTC")
		  )

return(st)
}
#' Draws sites of residency and adds a convex hull
#'
#' Draw a map (from the \code{R} Package \code{maps}) showing the defined
#' stationary sites
#'
#'
#' @param coord a \code{SpatialPoints} or \code{matrix} object, containing x
#' and y coordinates (in that order).
#' @param site a \code{numerical vector} assigning each row to a particular
#' period. Stationary periods in numerical order and values >0,
#' migration/movement periods 0.
#' @param points \code{logical}; if \code{TRUE}, the points of each site will
#' also be plottet.
#' @param map.range some possibilities to choose defined areas ("World
#' (default)", "EuroAfrica","America","AustralAsia").
#' @param xlim two element numeric vector giving a range of longitudes,
#' expressed in degrees, to which the map is restricted. Longitude is measured
#' in degrees east of Greenwich, so that, in particular, locations in
#' Switzerland have positive longitude. If map.range is defined this argument
#' will not be considered.
#' @param ylim two element numeric vector giving a range of latitudes,
#' expressed in degrees, to which the map is restricted. Latitude is measured
#' in degrees north of the equator, so that, in particular, locations in
#' Switzerland have positive latitude. If map.range is defined this argument
#' will not be considered.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param lwd The line width, a positive number, defaulting to 1.
#' @param lty The line type. Line types can either be specified as an integer
#' (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash,
#' 6=twodash) or as one of the character strings "blank", "solid", "dashed",
#' "dotted", "dotdash", "longdash", or "twodash", where "blank" uses "invisible
#' lines" (i.e., does not draw them).
#' @param pch Either an integer specifying a symbol or a single character to be
#' used as the default in plotting points. See \code{\link{points}} for
#' possible values and their interpretation. Note that only integers and
#' single-character strings can be set as a graphics parameter (and not NA nor
#' NULL).
#' @param cex A numerical value giving the amount by which plotting symbols
#' should be magnified relative to the default.
#' @param col a vector of colors with the same length as the number of defined
#' sites (if Default a predefined color ramp will be used).
#' @param main map title.
#' @param add \code{logical}; if \code{TRUE}, positions will be added to an
#' existing plot.
#' @author Simeon Lisovski
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' coord <- coord(tFirst,tSecond,type,-6)
#' filter <- distanceFilter(tFirst,tSecond,type,distance=30)
#' site <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)$site
#' siteMap(coord[filter,],site[filter],xlim=c(-20,20),ylim=c(0,60),lwd=2,pch=20,cex=0.5,main="hoopoe2")
#'
#' @export siteMap
#' @importFrom maps map
siteMap <-
function(coord,site,points=TRUE,map.range=c("EuroAfrica","AustralAsia","America","World"),xlim=NULL,ylim=NULL,xlab="Longitude",ylab="Latitude",lwd=1,lty=1,pch=1,cex=1,col="black",main=NULL,add=FALSE) {

nr.sites <- length(levels(as.factor(site[site!=0])))
site <- as.factor(site)


	if(sum(map.range==c("EuroAfrica","AustralAsia","America","World"))==4) {
		range <- c(-180,180,-75,90)
		} else
		{
			if(map.range %in% c("EuroAfrica","AustralAsia","America","World")) {
				if(sum(map.range=="EuroAfrica") ==1)    range <- c(-24,55,-55,70)
				if(sum(map.range=="AustralAsia") ==1)   range <- c(60,190,-55,78)
				if(sum(map.range=="America") ==1)       range <- c(-170,-20,-65,78)
				if(sum(map.range=="World") ==1)      	range <- c(-180,180,-75,90)
				} else {
				range <- c(-180,180,-75,90)
				}
		}

if(length(xlim)==0) range[1:2] <- range[1:2] else range[1:2] <- xlim
if(length(ylim)==0) range[3:4] <- range[3:4] else range[3:4] <- ylim

if(length(main)==0) main <- "" else main <- main


# Colors
colors <- c("firebrick1","forestgreen","orange","cornflowerblue","darkred","darkolivegreen3","darkorchid","black","chartreuse1","darkgrey","cyan2")


if(!add) {
	par(oma=c(5,3,0.5,0.5))
	map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),fill=T,lwd=0.01,col=c("grey90"),add=F,mar=c(rep(0.5,4)))
	map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),interior=TRUE,col=c("darkgrey"),add=TRUE)
	mtext(xlab,side=1,line=2.2,font=3)
	mtext(ylab,side=2,line=2.5,font=3)
	map.axes()

	mtext(main,line=0.6,cex=1.2)
	}

if(length(col)==length(levels(site))) col <- col else col <- colors[1:length(levels(site))]

sites <- length(levels(site[site!=0]))

if(points){
	for(i in 1:sites){
		points(coord[site==i,],cex=cex,pch=pch,col=col[i])
	}
	}

for(j in 1:sites){
      X <- na.omit(coord[site==j,])

      hpts <- chull(na.omit(coord[site==j,]))
      hpts <- c(hpts,hpts[1])
      lines(coord[hpts,],lty=lty,lwd=lwd,col=col[j])
    }

legend("bottomright",c(letters[1:nr.sites]),pch=pch, col=colors[1:nr.sites])

}

#' Write a file which plots a trip in Google Earth
#'
#' This function creates a .kml file from light intensity measurements over
#' time that can ve viewed as a trip in Google Earth.
#'
#'
#' @param file A character expression giving the whole path and the name of the
#' resulting output file including the .kml extension.
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
#' respectively
#' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
#' twilight"). Either a single value, a \code{vector} with the same length as
#' \code{tFirst}.
#' @param col.scheme the color scheme used for the points. Possible color
#' schemes are: \code{\link{rainbow}}, \code{\link{heat.colors}},
#' \code{\link{topo.colors}}, \code{\link{terrain.colors}}.
#' @param point.alpha a \code{numerical value} indicating the transparency of
#' the point colors on a scale from 0 (transparent) to 1 (opaque).
#' @param cex \code{numerical value} for the size of the points.
#' @param line.col An character expression (any of \code{\link{colors}} or
#' hexadecimal notation), or numeric indicating the color of the line
#' connecting the point locations.
#' @return This function returns no data. It creates a .kml file in the in the
#' defined path.
#' @author Simeon Lisovski and Michael U. Kemp
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' filter <- distanceFilter(tFirst,tSecond,type,distance=30)
#' trip2kml("trip.kml", tFirst[filter], tSecond[filter], type[filter],
#' 		degElevation=-6, col.scheme="heat.colors", cex=0.7,
#' 		line.col="goldenrod")
#'
#' @export trip2kml
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
#' @param coord a \code{SpatialPoints} or \code{matrix} object, containing x
#' and y coordinates (in that order).
#' @param equinox logical; if \code{TRUE}, the equinox period(s) is shown as a
#' broken blue line.
#' @param map.range some possibilities to choose defined areas (default:
#' "World").
#' @param xlim two element numeric vector giving a range of longitudes,
#' expressed in degrees, to which the map is restricted. Longitude is measured
#' in degrees east of Greenwich, so that, in particular, locations in
#' Switzerland have positive longitude. If map.range is defined this argument
#' will not be considered.
#' @param ylim two element numeric vector giving a range of latitudes,
#' expressed in degrees, to which the map is restricted. Latitude is measured
#' in degrees north of the equator, so that, in particular, locations in
#' Switzerland have positive latitude. If map.range is defined this argument
#' will not be considered.
#' @param xlab a title for the x axis.
#' @param ylab a title for the y axis.
#' @param col colour code or name for the trip line (default equals red).
#' @param main map title
#' @param legend \code{logical}; if \code{TRUE}, a legend will be added to the
#' plot.
#' @param add \code{logical}; if \code{TRUE}, positions will be added to an
#' existing plot.
#' @author Simeon Lisovski
#' @examples
#'
#' data(hoopoe2)
#' attach(hoopoe2)
#' coord <- coord(tFirst,tSecond,type,-6)
#' tripMap(coord,xlim=c(-20,20),ylim=c(0,60),main="hoopoe2")
#'
#' @export tripMap
tripMap <-
function(coord,equinox=TRUE,map.range=c("EuroAfrica","AustralAsia","America","World"),xlim=NULL,ylim=NULL,xlab="Longitude",ylab="Latitude",col="tomato2",main=NULL,legend=TRUE,add=FALSE) {


	if(sum(map.range==c("EuroAfrica","AustralAsia","America","World"))==4) {
		range <- c(-180,180,-75,90)
		} else
		{
			if(map.range %in% c("EuroAfrica","AustralAsia","America","World")) {
				if(sum(map.range=="EuroAfrica") ==1)    range <- c(-24,55,-55,70)
				if(sum(map.range=="AustralAsia") ==1)   range <- c(60,190,-55,78)
				if(sum(map.range=="America") ==1)       range <- c(-170,-20,-65,78)
				if(sum(map.range=="World") ==1)      	range <- c(-180,180,-75,90)
				} else {
				range <- c(-180,180,-75,90)
				}
		}

if(length(xlim)==0) range[1:2] <- range[1:2] else range[1:2] <- xlim
if(length(ylim)==0) range[3:4] <- range[3:4] else range[3:4] <- ylim

if(length(main)==0) main <- "" else main <- main


if(!add) {
	par(oma=c(5,3,0,0))
	map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),fill=T,lwd=0.01,col=c("grey90"),add=F,mar=c(rep(0.5,4)))
	map(xlim=c(range[1],range[2]),ylim=c(range[3],range[4]),interior=TRUE,col=c("darkgrey"),add=TRUE)
	mtext(xlab,side=1,line=2.2,font=3)
	mtext(ylab,side=2,line=2.5,font=3)
	map.axes()

	mtext(main,line=0.6,cex=1.2)
	}

	points(coord,pch=3,cex=0.7)
	lines(coord,lwd=0.5,col=col)

if(equinox){
	nrow <- 1
	repeat{
		while(is.na(coord[nrow,2])==FALSE) {
			nrow <- nrow + 1
			if(nrow==nrow(coord)) break
			}
			if(nrow==nrow(coord)) break
			start   <- nrow-1
		while(is.na(coord[nrow,2])) {
			nrow <- nrow + 1
			if(nrow==nrow(coord)) break
			}
			if(nrow==nrow(coord)) break
			end    <- nrow

		lines(c(coord[start,1],coord[end,1]),c(coord[start,2],coord[end,2]),col="blue",lwd=3,lty=1)
		}
if(legend) legend("bottomright",lty=c(0,1,1),pch=c(3,-1,-1),lwd=c(1,0.5,3),col=c("black",col,"blue"),c("Positions","Trip","Equinox"),bty="n",bg="grey90",border="grey90",cex=0.8)
} else {
	if(legend) legend("bottomright",lty=c(0,1),pch=c(3,-1),lwd=c(1,0.5),col=c("black",col,"blue"),c("Positions","Trip"),bty="n",bg="grey90",border="grey90",cex=0.8)
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
twilightCalc <- function(datetime, light, LightThreshold=TRUE, preSelection=TRUE, maxLight=NULL, ask=TRUE, nsee=500, allTwilights=FALSE)
{
  bas <- data.frame(datetime=as.POSIXct(as.character(datetime),"UTC"),light)

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

