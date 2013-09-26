# Function to summerize the individual migration pattern (Package GeoLight)
# ---------------------------------------------------------------------------
# Author: Simeon Lisovski, November 2011
# ---------------------------------------------------------------------------



#' Summary of migration/movement pattern
#' 
#' Function for making a data frame summarising residency and movement pattern.
#' 
#' 
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param site a \code{vector}, indicating the residency period of a particular
#' day (see output: \code{\link{changeLight}})
#' @return A \code{data.frame} with end and start date (yyyy-mm-dd hh:mm, UTC)
#' for each stationary period.
#' @author Simeon Lisovski
#' @examples
#' 
#' data(hoopoe2)
#' attach(hoopoe2)
#' residency <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)
#' schedule(tFirst,tSecond,residency$site)
#' 
#' @export schedule
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
