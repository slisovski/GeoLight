\name{HillEkstromCalib}
\alias{HillEkstromCalib}
\title{Hill-Ekstrom calibration}
\description{Hill-Ekstrom calibration for one or multiple stationary periods.}
\usage{
HillEkstromCalib(tFirst,tSecond,type,site, 
				start.angle=-6,distanceFilter=FALSE,
				distance,plot=TRUE)
}

\arguments{
  \item{tFirst}{date and time of sunrise/sunset (e.g. 2008-12-01 08:30)}
  \item{tSecond}{date and time of sunrise/sunset (e.g. 2008-12-01 17:30)}
  \item{type}{either 1 or 2, defining \code{tFirst} as sunrise or sunset respectively}
  \item{site}{a \code{numerical vector} assigning each row to a particular period. Stationary periods in numerical order and values >0, migration/movement periods 0}
  \item{start.angle}{a single sun elevation angle. The combined process of checking for minimal variance in resulting latitude, which is the initial value for the sun elevation angle in the iterative process of identifying the latitudes with the least variance}
  \item{distanceFilter}{logical, if TRUE the \code{\link{distanceFilter}} will be used to filter unrealistic positions}
  \item{distance}{if \code{distanceFilter} is set \code{TRUE} a threshold distance in km has to be set (see: \code{\link{distanceFilter}})}
  \item{plot}{logical, if TRUE the function will give a plot with all relevant information}
}
\value{A \code{vector} of sun elevation angles corresponding to the Hill-Ekstrom calibration for each defined period.}

\details{The \emph{Hill-Ekstrom calibration} has been suggested by Hill & Braun (2001) and Ekstrom (2004), and allows for calibrating data during stationary periods at unknown latitudinal positions. The Hill-Ekstrom calibration bases on an increasing error range in latitudes with an increasing mismatch between light level threshold and the used sun angle. This error is strongly amplified with proximity to the equinox times due to decreasing slope of day length variation with latitude. Furthermore, the sign of the error switches at the equinox, i.e. latitude is overestimated before the equinox and underestimated after the equinox (or vice versa depending on autumnal/vernal equinox, hemisphere, and sign of the mismatch between light level threshold and sun angle). When calculating the positions of a stationary period, the variance in latitude is minimal if the sun elevation angle fits to the defined light level threshold. Moreover, the accuracy of positions increases with decreasing variance in latitudes. \bold{However, the method is only applicable for stationary periods and under stable shading intensities}. The plot produced by the function may help to judge visually if the calculated sun elevation angles are realistic (e.g. site 2 in the example below) or not (e.g. site 3 in the example below. }
\references{
Ekstrom, P.A. (2004) An advance in geolocation by light. \emph{Memoirs of the National Institute of Polar Research}, Special Issue, \bold{58}, 210-226.

Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step: Latitude. In: \emph{Electronic Tagging and Tracking in Marine Fisheries} (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The Netherlands.

Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and precision affected by environmental factors. \emph{Methods in Ecology and Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.}

\author{Simeon Lisovski}
\examples{
data(hoopoe2)
attach(hoopoe2)
residency <- changeLight(tFirst,tSecond,type,rise.prob=0.1,set.prob=0.1,plot=FALSE,summary=FALSE)
HillEkstromCalib(tFirst,tSecond,type,residency$site,-6)
}