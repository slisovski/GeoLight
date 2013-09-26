\name{lightFilter}
\alias{lightFilter}
\title{Filter to remove noise in light intensity measurements during the night}
\description{The filter identifies and removes light intensities oczillating around the baseline or few light intensities resulting in a short light peak during the night. Such noise during the night will increase the calculated twilight events using the function \code{\link{twilightCalc}} and therewith the manual work to remove these false twilight events.}
\usage{lightFilter(light, baseline=NULL, iter=2)}
\arguments{
   \item{light}{\code{numerical} value of the light intensity (usually arbitrary units).}
  \item{baseline}{the light intensity baseline (no light). If \code{Default}, it will be calculated as the most frequent value below the mean light intensities.}
  \item{iter}{a \code{numerical} value, specifying how many iterations should be computed (see details).}
}
\value{numerical \code{vector} with the new light levels. Same length as the initial light vector.}
\details{The filter searches for light levels above the baseline and compares the prior and posterior levels. If these values are below the threshold the particular light level will be reduced to the baseline. A few (usually two) iterations might be enough to remove most noise during the night (however, not if such noise occurs at the begining or at the end were not enough prior or posterior values are available).}
\author{Simeon Lisovski}
\examples{
night <- rep(0,50); night[runif(4,0,50)] <- 10; night[runif(4,0,50)] <- -5
nightday <- c(night,rep(30,50))
plot(nightday,type="l",ylim=c(-5,30),ylab="light level",xlab="time (time)")
light2 <- lightFilter(nightday, baseline=0, iter=4)
lines(light2,col="red")
legend("bottomright",c("before","after"),lty=c(1,1),col=c("black","red"),bty="n")
}