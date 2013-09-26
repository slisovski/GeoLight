\name{schedule}
\alias{schedule}
\title{Summary of migration/movement pattern}

\description{Function for making a data frame summarising residency and movement pattern.}

\usage{schedule(tFirst,tSecond,site)}

\arguments{
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{site}{a \code{vector}, indicating the residency period of a particular day (see output: \code{\link{changeLight}})}
}
\value{A \code{data.frame} with end and start date (yyyy-mm-dd hh:mm, UTC) for each stationary period.}

\author{
Simeon Lisovski
}

\examples{
data(hoopoe2)
attach(hoopoe2)
residency <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)
schedule(tFirst,tSecond,residency$site)
}
