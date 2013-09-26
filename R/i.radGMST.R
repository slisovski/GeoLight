i.radGMST <-
function(jD,jD0,jC,jC0) {

#--------------------------------------------------------------------------------------------------------
# jD:  Julian Date with Hour and Minute
# jD0: Julan Date at t 0 UT
# jC:  Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
# jC0: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00) at t 0 UT
#--------------------------------------------------------------------------------------------------------

options(digits=10)


		UT  <- 86400 * (jD-jD0)
		st0 <- 24110.54841 + 8640184.812866*jC0 + 1.0027379093*UT + (0.093104 - 0.0000062*jC0)*jC0*jC0
		gmst<- (((2*pi)/86400)*(st0%%86400))

return(gmst)

}

