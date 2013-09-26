# Summery of the changeLight Function
# ---------------------------------------------------------------------------
# Author: Simeon Lisovski, Mai 2012
# ---------------------------------------------------------------------------

i.sum.Cl <- function(object) {
	
	if(sum(names(object)==c("riseProb","setProb","rise.prob","set.prob","site","migTable"))==6){
		cat("\n")
		cat("Probability threshold(s):")
		cat(rep("\n",2))
		if(!is.na(object$rise.prob)) cat(paste("	Sunrise: ",object$rise.prob))
		if(!is.na(object$set.prob)) cat(paste("	Sunset: ",object$set.prob))
		cat(rep("\n",3))
		cat("Migration schedule table:")
		cat(rep("\n",2))
		
		print(object$migTable,quote=FALSE)
	} else {
		cat("Error: List must be the output list of the changeLight function.")
	} 
}