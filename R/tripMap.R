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


if(add==FALSE) {	
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

if(equinox==TRUE){
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
if(legend==TRUE) legend("bottomright",lty=c(0,1,1),pch=c(3,-1,-1),lwd=c(1,0.5,3),col=c("black",col,"blue"),c("Positions","Trip","Equinox"),bty="n",bg="grey90",border="grey90",cex=0.8)
} else {
	if(legend==TRUE) legend("bottomright",lty=c(0,1),pch=c(3,-1),lwd=c(1,0.5),col=c("black",col,"blue"),c("Positions","Trip"),bty="n",bg="grey90",border="grey90",cex=0.8)
	}	

}

