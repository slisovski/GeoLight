i.preSelection <- function(datetime, light, LightThreshold){
	
	dt <- cut(datetime,"1 hour")
	st <- as.POSIXct(levels(dt),"UTC")
	
	raw <- data.frame(datetime=dt,light=light)
	
	h  <- tapply(light,dt,max)
	df1 <- data.frame(datetime=st+(30*60),light=as.numeric(h))
	
    smooth <- i.twilightEvents(df1[,1], df1[,2], LightThreshold)
    		  smooth <- data.frame(id=1:nrow(smooth),smooth)
	raw    <- i.twilightEvents(datetime, light, LightThreshold)
			  raw <- data.frame(id=1:nrow(raw),raw)
			
  ind2 <- rep(NA,nrow(smooth))
  for(i in 1:nrow(smooth)){
  	tmp <- subset(raw,datetime>=(smooth[i,2]-(90*60)) & datetime<=(smooth[i,2]+(90*60)))
  	
  	if(smooth[i,3]==1) ind3 <- tmp$id[which.min(tmp[,2])]
  	if(smooth[i,3]==2) ind3 <- tmp$id[which.max(tmp[,2])]
  ind2[i] <- ind3
  }


res <- data.frame(raw,mod=1)
res$mod[ind2] <- 0

return(res)
}