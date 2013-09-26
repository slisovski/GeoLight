i.radObliquity <-
function(jC) {

#--------------------------------------------------------------------------------------------------------
# jC: Number of julian centuries from the julianian epoch J2000 (2000-01-01 12:00)
#--------------------------------------------------------------------------------------------------------

options(digits=10)


degObliquity <- 23.43929111 - (46.8150 + (0.00059 - 0.001813 *jC)*jC) *jC/3600
radObliquity <- i.rad(degObliquity)

return(radObliquity)		

}

