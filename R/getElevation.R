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
 
 if(plot==TRUE)
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