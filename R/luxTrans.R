luxTrans <-
function(file="/path/file.lux") {
	
lux1 <- read.table(file,sep="\t",skip=21,col.names=c("datetime","time")) # read file
lux <- data.frame(datetime=as.POSIXct(strptime(lux1[,1],format="%d/%m/%Y %H:%M:%S",tz="UTC")),light=lux1[,2])

return(lux)
}

