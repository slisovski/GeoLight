i.get.outliers<-function(residuals, k=3) {
	x <- residuals
	# x is a vector of residuals
	# k is a measure of how many interquartile ranges to take before saying that point is an outlier
	# it looks like 3 is a good preset for k
	QR<-quantile(x, probs = c(0.25, 0.75))
	IQR<-QR[2]-QR[1]
	Lower.band<-QR[1]-(k*IQR)
	Upper.Band<-QR[2]+(k*IQR)
	delete<-which(x<Lower.band |  x>Upper.Band)
	return(as.vector(delete))
}

