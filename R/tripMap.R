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

