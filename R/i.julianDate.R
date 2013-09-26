i.julianDate <-
function(year,month,day,hour,min) {

#--------------------------------------------------------------------------------------------------------
# Year:		Year as numeric e.g. 2010
# Month:	Month as numeric e.g.   1
# Day:		Day as numeric e.g.		1
# Hour:		Hour as numeric e.g.   12
# Min:		Minunte as numeric e.g. 0
#--------------------------------------------------------------------------------------------------------

	options(digits=15)

			fracOfDay	<- hour/24 + min/1440
		
			# Julian date (JD)
			# ------------------------------------
      
      		index1 <- month <= 2
			if(sum(index1) > 0)
				{
					year[index1]  <- year[index1] -1
					month[index1] <- month[index1] +12
				}
			
     		index2 <- (year*10000)+(month*100)+day <= 15821004
     
      		JD <- numeric(length(year))	
			if(sum(index2)>0)
				{
					JD[index2] <- floor(365.25*(year[index2]+4716)) + floor(30.6001*(month[index2]+1)) + day[index2] + fracOfDay[index2] - 1524.5
				} 
		  	index3 <- year*10000+month*100+day >= 15821015
      		if (sum(index3)>0)
				{
					a <- floor(year/100)
					b <- 2 - a + floor(a/4)
					
					JD[index3] <- floor(365.25*(year[index3]+4716)) + floor(30.6001*(month[index3]+1)) + day[index3] + fracOfDay[index3] + b[index3] - 1524.5
				} 
        	JD[!index2&!index3] <- 1 
		
return(JD)
				
}

