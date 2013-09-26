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
	

cat(paste("Note: ",length(index[index==FALSE])," of ",length(index[coord[,2]!=999])," positions were filtered (",floor((length(index[index==FALSE])*100)/length(index[coord[,2]!=999]))," %)",sep=""))
return(index)
}
