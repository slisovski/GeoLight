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
  
  if(preSelection==FALSE) out$mod <- 0
  
  if(ask==TRUE) 	
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
 
  if(allTwilights==TRUE) {
  	results$consecTwilights <- opt
  	return(results)
  } else {
  return (opt)}
}