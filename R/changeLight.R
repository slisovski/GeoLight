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
	CPs1 <- binseg.mean.cusum(rise, Q=length(rise)/2, pen=0.001)
	CPs2 <- binseg.mean.cusum(set, Q=length(set)/2, pen=0.001) 

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
	
if(plot==T){
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
	if(comp$filt[1]==FALSE) site[1] <- 1
	for(i in 2:length(date)){
		
	if(comp$filt[i]==FALSE & site[i-1]!=0) site[i] <- as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])
	if(comp$filt[i]==FALSE & site[i-1]==0) site[i] <- as.numeric(levels(as.factor(site))[length(levels(as.factor(site)))])+1
	if(comp$filt[i]==TRUE & abs(as.numeric(
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
	

if(plot==T) {
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


if(summary==TRUE){i.sum.Cl(out)}

return(out)

} # end function