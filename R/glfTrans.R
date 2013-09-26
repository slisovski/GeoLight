glfTrans <-
function(file="/path/file.glf") {

	
glf1 <- read.table(file,sep="\t",skip=6,col.names=c("datetime","light","1","2","3")) # read file
 
# Date transformation
 	year   <- as.numeric(substring(glf1$datetime,7,10))
 	month  <- as.numeric(substring(glf1$datetime,4,5))
 	day    <- as.numeric(substring(glf1$datetime,1,2))
 	hour   <- as.numeric(substring(glf1$datetime,12,13))
 	min    <- as.numeric(substring(glf1$datetime,15,16))
 	gmt.date <- paste(year,"-", month,"-",day," ",hour,":",min,":",0,sep="")
 	gmt.date <- as.POSIXct(strptime(gmt.date, "%Y-%m-%d %H:%M:%S"), "UTC")

glf <- data.frame(datetime=gmt.date,light=glf1$light)

return(glf)
}

