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


if(add==FALSE) {	
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

if(points==TRUE){
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

