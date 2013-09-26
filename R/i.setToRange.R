i.setToRange <-
function(Start,Stop,Angle) {

#-------------------------------------------------------------------------------------------------------------------
# Start:	 Minimal value of the range in degrees
# Stop: 	 Maximal value of the range in degrees
# Angle:	 The angle that should be fit to the range
#-------------------------------------------------------------------------------------------------------------------

options(digits=15)
		
		angle <- Angle
		range <- Stop - Start
				
			
			
			index1 <- angle >= Stop
			if (sum(index1)>0) angle[index1] <- angle[index1] - (floor((angle[index1]-Stop)/range)+1)*range
			
			index2 <- angle < Start
			if(sum(index2)>0) angle[index2] <- angle[index2]  + ceiling(abs(angle[index2] -Start)/range)*range
				
return(angle)
			
}

