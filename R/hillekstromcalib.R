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

