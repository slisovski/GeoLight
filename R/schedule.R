# Function to summerize the individual migration pattern (Package GeoLight)
# ---------------------------------------------------------------------------
# Author: Simeon Lisovski, November 2011
# ---------------------------------------------------------------------------

schedule <- function(tFirst,tSecond,site){
	
	tFirst <- as.POSIXct(as.character(tFirst),"UTC")
	tSecond <- as.POSIXct(as.character(tSecond),"UTC")
	
	
	midnoon <- tFirst + (tSecond - tFirst)/2
	diff1 <- c(site,site[length(site)]) != c(0,site)
	a <- which(diff1)[c(T,F)]
	b <- which(diff1)[c(F,T)]

	rows <- sort(c(
					(if(a[1]==1) c(a[1],a[-1]-1) else a-1)
					,b,length(midnoon)))
	
	midnoon1 <- midnoon[rows]
	
	st <- data.frame(site=c(letters[1:(length(midnoon1)/2)]),
		  start=as.POSIXct(midnoon1[c(T,F)],tz="UTC"),
		  end=as.POSIXct(midnoon1[c(F,T)],tz="UTC")
		  )

return(st)
}