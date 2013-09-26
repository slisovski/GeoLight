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