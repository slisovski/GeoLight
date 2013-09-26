i.twilightEvents <- function(datetime, light, LightThreshold){
   
   df <- data.frame(datetime, light)
   
   ind1 <- which((df$light[-nrow(df)] < LightThreshold & df$light[-1] > LightThreshold) | 
    			 (df$light[-nrow(df)] > LightThreshold & df$light[-1] < LightThreshold) | 
  				  df$light[-nrow(df)] == LightThreshold)

   bas1 <- cbind(df[ind1,],df[ind1+1,])
  		  bas1 <- bas1[bas1[,2]!=bas1[,4],]
  
  x1 <- as.numeric(unclass(bas1[,1])); x2 <- as.numeric(unclass(bas1[,3]))
  y1 <- bas1[,2]; y2 <- bas1[,4]
  m <- (y2-y1)/(x2-x1)
  b <- y2-(m*x2)
  
  xnew <- (LightThreshold - b)/m
  type <- ifelse(bas1[,2]<bas1[,4],1,2)
  res  <- data.frame(datetime=as.POSIXct(xnew, origin="1970-01-01", tz="UTC"),type)
  
return(res)
	
}