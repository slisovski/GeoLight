i.radDeclination <-
function(radEclipticalLength,radObliquity) {

#-------------------------------------------------------------------------------------------------------------------
# RadEclipticLength: The angle between an object's rotational axis, and a line perpendicular to its orbital plane. 
# EadObliquity:
#-------------------------------------------------------------------------------------------------------------------

options(digits=10)

	dec <- asin(sin(radObliquity)*sin(radEclipticalLength))
			  	
return(dec)

}

