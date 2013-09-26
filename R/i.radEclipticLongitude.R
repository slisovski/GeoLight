i.radEclipticLongitude <-
function(jC) {

#-------------------------------------------------------------------------------------------------------------------
# jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00
#-------------------------------------------------------------------------------------------------------------------



	options(digits=10)

  	radMeanAnomaly <- 2*pi*i.frac(0.993133 + 99.997361*jC)
	EclipticLon    <- 2*pi*i.frac(0.7859452 + radMeanAnomaly/(2*pi) + (6893*sin(radMeanAnomaly) + 72*sin(2*radMeanAnomaly) + 6191.2*jC) / 1296000)

return(EclipticLon)

}

