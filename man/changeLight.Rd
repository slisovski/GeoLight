\name{changeLight}
\alias{changeLight}
\title{Residency analysis using a changepoint model}

\description{Function to discriminate between periods of residency and movement based on consecutive sunrise and sunset data. The calculation is based on a changepoint model (\bold{\pkg{R}} Package \code{\link{changepoint}}: \code{\link{binseg.mean.cusum}}) to find multiple changepoints within the data.}

\usage{changeLight(tFirst,tSecond,type,quantile=0.6,rise.prob=NA,set.prob=NA,days=5,plot=TRUE,summary=TRUE)}

\arguments{
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively}
  \item{quantile}{probability threshold for stationary site selection. Higher values (above the defined quantile of all probabilities) will be considered as changes in the behavior. If specified, \code{rise.prob} and/or \code{set.prob} will not be considered.}
  \item{rise.prob}{the probability threshold for \bold{sunrise}: greater or equal values indicates changes in the behaviour of the individual. If, \code{NA} sunrise will not be considered.}
  \item{set.prob}{the probability threshold for \bold{sunset}: higher and equal values indicates changes in the behaviour of the individual. If, \code{NA} sunrise will not be considered.}
  \item{days}{a threshold for the length of stationary period. Periods smaller than "days" will not be considered as a residency period}
  \item{plot}{logical, if \code{TRUE} a plot will be produced}
  \item{summary}{logical, if \code{TRUE} a summary of the results will be printed}
}
\details{The \code{binseg.mean.cusum} from the code{R} Package \code{changepoint} is a function to find a multiple changes in mean for data where no assumption is made on their distribution. The value returned is the result of finding the optimal location of up to Q changepoints (in this case as many as possible) using the cumulative sums test statistic.}

\value{A \code{list} with probabilities for \emph{sunrise} and \emph{sunset} the user settings of the probabilities and the resulting stationary periods given as a \code{vector}, with the residency sites as positiv numbers in ascending order (0 indicate movement/migration).}

\note{The sunrise and/or sunset times shown in the graph (if \code{plot=TRUE}) represent hours of the day. However if one or both of the twilight events cross midnight during the recording period the values will be transformed to avoid discontinuity.}

\seealso{\code{\link{changepoint}}, \code{\link{binseg.mean.cusum}}}

\references{
Taylor, Wayne A. (2000) Change-Point Analysis: A Powerful New Tool For Detecting Changes.

M. Csorgo, L. Horvath (1997) Limit Theorems in Change-Point Analysis. \emph{Wiley}.
	
Chen, J. and Gupta, A. K. (2000) Parametric statistical change point analysis. \emph{Birkhauser}.
}
\author{Simeon Lisovski & Tamara Emmenegger}
\examples{
data(hoopoe2)
attach(hoopoe2)
residency <- changeLight(tFirst,tSecond,type,quantile=0.9)
}