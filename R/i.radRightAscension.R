i.radRightAscension <-
function(RadEclipticalLength,RadObliquity) {

#-------------------------------------------------------------------------------------------------------------------
# RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane. 
# EadObliquity:
#-------------------------------------------------------------------------------------------------------------------

options(digits=10)

	index1 <- (cos(RadEclipticalLength) < 0)
	
	res <- numeric(length(RadEclipticalLength))
	if (sum(index1)>0) 
		{ 
			res[index1] <- (atan((cos(RadObliquity[index1])*sin(RadEclipticalLength[index1]))/cos(RadEclipticalLength[index1])) + pi) 
		}
	
	index2 <- (cos(RadEclipticalLength) >= 0)
	if (sum(index2)>0)
		{ 
			res[index2] <- (atan((cos(RadObliquity[index2])*sin(RadEclipticalLength[index2]))/cos(RadEclipticalLength[index2]))) 
		}
		

return(res)

}

