\name{calib*}
\alias{calib1}
\alias{calib2}
\alias{calib}
\title{
Example data for calibration: Light intensities and twilight events
}
\description{
Light intensity measurements over time (calib1) recorded at the rooftop of the Swiss Ornithological Institute (Lon: 8.0, Lat: 47.01). Defined twilight events from calib1 (calib2). These data serve as an example for calculating the sun elevation angle of an additional data set, which is subsequently used to calibrate the focal dataset.
}
\usage{
data(calib1)

data(calib2)
}
\references{Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and precision affected by environmental factors. \emph{Methods in Ecology and Evolution}, DOI: 10.1111/j.2041-210X.2012.00185.x.
}
\examples{
data(calib2)
attach(calib2)
getElevation(tFirst,tSecond,type,c(8,47.01))
}