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

if(plot==TRUE){	
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
all <- rbind(subset(dusk,filter==TRUE),subset(dawn,filter==TRUE))

filter <- rep(FALSE,length(tFirst))
	filter[tFirst%in%all$datetime | tSecond%in%all$datetime] <- TRUE
	
return(!filter)
}