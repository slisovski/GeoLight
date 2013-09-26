lightFilter <- function(light, baseline=NULL, iter=2){
	
	r <- as.data.frame(table(light))
    	 r[,1] <- as.character(r[,1])
    nr <- as.numeric(which.max(r$Freq[as.numeric(r[,1])<mean(light)]))
    LightThreshold <- ifelse(is.null(baseline),as.numeric(r[nr,1]),baseline)
    
    light[light<LightThreshold] <- LightThreshold
    
	
	index <- which(light<mean(light) & light!= LightThreshold)
	
	rep   <- rep(FALSE,length(light))
	
	for(i in 1:iter){
		
	for(i in index[index>5 & index<(length(light)-5)]){
		
	back=FALSE
	if(any(light[seq(i-5,i)]==LightThreshold)) back <- TRUE	
	forw=FALSE	
	if(any(light[seq(i,i+5)]==LightThreshold)) forw <- TRUE
	
	if(back==TRUE & forw==TRUE) rep[i] <- TRUE
	
	}

	light[rep] <- LightThreshold
	
	}
	
return(light)	
	
}