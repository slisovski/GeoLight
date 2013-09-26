gleTrans <- function(file) {


gle1 <- read.table(file,sep="\t",skip=16,col.names=c("date","light","1","2","3","4","5","6","7","8","type")) # read file
 	    gle1 <- subset(gle1,gle1$type>0,select=c("date","type"))
 
# Date transformation
 	year   <- as.numeric(substring(gle1$date,7,10))
 	month  <- as.numeric(substring(gle1$date,4,5))
 	day    <- as.numeric(substring(gle1$date,1,2))
 	hour   <- as.numeric(substring(gle1$date,12,13))
 	min    <- as.numeric(substring(gle1$date,15,16))
 	gmt.date <- paste(year,"-", month,"-",day," ",hour,":",min,":",0,sep="")
 	gmt.date <- as.POSIXct(strptime(gmt.date, "%Y-%m-%d %H:%M:%S"), "UTC")

gle <- data.frame(date=gmt.date,type=gle1$type)


opt <- data.frame(tFirst=as.POSIXct("1900-01-01 01:01","UTC"),tSecond=as.POSIXct("1900-01-01 01:01","UTC"),type=0)

row <- 1
for (i in 1:(length(gmt.date)-1))
	{
	  if (abs(as.numeric(difftime(gle$date[i],gle$date[i+1],units="hours")))< 18 & gle$date[i] != gle$date[i+1])
	  	{
	  		opt[row,1] <- gle$date[i]
	  		opt[row,2] <- gle$date[i+1]
			opt[row,3] <- gle$type[i]
			
			row <- row+1
		}
	}

return(opt)
}

