\name{GeoLight-package}
\alias{GeoLight}
\title{The GeoLight Package}

\description{
This is a summary of all features of \bold{\code{GeoLight}}, a \code{R}-package for analyzing light based geolocator data.
}

\section{Details}{
\bold{\code{GeoLight}} is a package to derive geographical positions from daily light intensity pattern. Positioning and calibration methods are based on the threshold-method (Ekstrom 2004, Lisovski \emph{et al.} 2012). A changepoint model from the \code{R} package \code{changepoint} is implemented to distinguish between periods of residency and movement based on the sunrise and sunset times. Mapping functions are implemented using the \code{R} package \code{maps}.
}

\section{Getting Started}{
We refrain from giving detailed background on the (several steps of) analysis of light-based geolocator data here but strongly recommend the key-publications below.
}

\section{Updates}{
We advise all users to update their installation of \bold{\code{GeoLight}} regularly.

Type \code{news(package="GeoLight")} to read news documentation about changes to the recent and all previous version of the package
}

\section{Important Notes}{
Most functions in \bold{\code{GeoLight}} require the same initial units and mostly the format and object type is mandatory:
 
    \tabular{rll}{
    \tab \code{tFirst} \tab yyyy-mm-dd hh:mm "UTC" (see: \code{\link{as.POSIXct}}, \link[=Sys.timezone]{time zones})\cr
	\tab \code{tSecond} \tab as \emph{tFirst} (e.g. 2008-12-01 17:30) \cr
	\tab \code{type} \tab either 1 or 2 depending on wheter \emph{tFirst} is sunrise (1) or sunset (2)\cr
	\tab \code{coord} \tab \code{SpatialPoints} or a \code{matrix}, containing x and y coordinates (in that order) \cr
	\tab \code{degElevation} \tab a \code{vector} or a single \code{value} of sun elevation angle(s) in degrees (e.g. -6)
       }
}


\section{FUNCTIONS AND DATASETS}{
In the following, we give a summary of the main functions and sample datasets in the \bold{\code{GeoLight}} package. Alternatively a list of all functions and datasets in alphabetical order is available \cr by typing \code{library(help=GeoLight)}. For further information on any of these functions, type \code{help(function name)}.
}


\section{CONTENTS:}{
\tabular{rll}{
\tab {I.} \tab Determination of sunset and sunrise \cr
\tab {II.} \tab Residency analysis \cr
\tab {III.} \tab Calibration \cr
\tab {IV.} \tab {Positioning}\cr
\tab {V.} \tab {Data visualisation}\cr
\tab {VI.} \tab {Examples}
}
}


\section{I. Determination of sunset and sunrise}{
\tabular{rll}{
\tab \code{\link{gleTrans}} \tab transformation of already defined twilight events* \cr
\tab \code{\link{glfTrans}} \tab transformation of light intensity measurements over time* \cr
\tab \code{\link{luxTrans}} \tab transformation of light intensity measurements over time** \cr
\tab \code{\link{lightFilter}} \tab filter to remove noise in light intensity measurements during the night \cr
\tab \code{\link{twilightCalc}} \tab definition of twilight events (\emph{sunrise, sunset}) from light intensity measurements \cr
}

* written for data recorded by geolocator devices from the \bold{Swiss Ornithological Institute}
** written for data recorded by geolocator devices from \bold{Migrate Technology Ltd}
}


\section{II. Residency Analysis}{
\tabular{rll}{
\tab \code{\link{changeLight}} \tab function to distinguish between residency and movement periods \cr
\tab \code{\link{schedule}} \tab function to produce a data frame summerizing the residency and movement pattern \cr
}
}


\section{III. Calibration}{

See Lisovski \emph{et al.} 2012 for all implemented calibration methods.

\tabular{rll}{
\tab \code{\link{getElevation}} \tab function to calculate the sun elevation angle for data with known position \cr
\tab \code{\link{HillEkstromCalib}} \tab \emph{Hill-Ekstrom calibration} for one or more defined stationary periods \cr
}
}


\section{IV. Positioning}{
\tabular{rll}{
\tab \code{\link{coord}} \tab main function to derive a \code{matrix} of spatial coordinates \cr
\tab \code{\link{distanceFilter}} \tab filter function to reduce unrealistic positions (not recommended, since the filtering ignore positioning error) \cr
\tab \code{\link{loessFilter}} \tab filter function to define outliers in sunrise and sunset times (defined twilight events) \cr
}
}


\section{V. Data visualisation}{
\tabular{rll}{
\tab \code{\link{tripMap}} \tab function to map the derived positions and combine the coordinates in time order\cr
\tab \code{\link{siteMap}} \tab function to show the results of the residency analysis on a map \cr
}
}


\section{IV. Examples}{
\tabular{rll}{
\tab \code{\link{calib1}} \tab data for calibration: light intensities \cr
\tab \code{\link{calib2}} \tab  data for calibration: Calculated twilight events (from \code{\link{calib1}} by \code{\link{twilightCalc}}) \cr
\tab \code{\link{hoopoe1}} \tab light intensity measurements over time recorded on a migratory bird \cr
\tab \code{\link{hoopoe2}} \tab sunrise and sunset times: From light intensity measurement (from \code{\link{hoopoe1}}) \cr
}
}


\section{\bold{\code{R}} Packages for Further Spatial Ananlyses}{
\code{spatstat} \cr
\code{adehabitat} \cr
\code{gstat} \cr
\code{trip} \cr
\code{tripEstimation} \cr
\code{move} \cr
...
}


\section{Acknowledgements}{
Steffen Hahn, Felix Liechti, Fraenzi Korner-Nievergelt, Andrea Koelzsch, Eldar Rakhimberdiev, Michael Sumner, Erich Baechler, Eli Bridge,  Andrew Parnell, Richard Inger
}

\section{Authors}{
Simeon Lisovski, Silke Bauer, Tamara Emmenegger

Maintainer: Simeon Lisovski: <\email{migration@vogelwarte.ch}>

Discussion Forum: \url{http://simeonlisovski.wordpress.com/geolight/}
}

\section{References}{
Ekstrom, P.A. (2004) An advance in geolocation by light. \emph{Memoirs of the National Institute of Polar Research}, Special Issue, \bold{58}, 210-226.

Fudickar, A.M., Wikelski, M., Partecke, J. (2011) Tracking migratory songbirds: accuracy of light-level loggers (geolocators) in forest habitats. \emph{Methods in Ecology and Evolution}, DOI: 10.1111/j.2041-210X.2011.00136.x.

Hill, C. & Braun, M.J. (2001) Geolocation by light level - the next step: Latitude. \emph{Electronic Tagging and Tracking in Marine Fisheries} (eds J.R. Sibert & J. Nielsen), pp. 315-330. Kluwer Academic Publishers, The Netherlands.

Hill, R.D. (1994) Theory of geolocation by light levels. \emph{Elephant Seals: Population Ecology, Behavior, and Physiology} (eds L. Boeuf, J. Burney & R.M. Laws), pp. 228-237. University of California Press, Berkeley.

Lisovski, S. and Hahn, S. (2012) GeoLight - processing and analysing light-based geolocator data in R. \emph{Methods in Ecology and Evolution}, doi: 10.1111/j.2041-210X.2012.00248.x

Lisovski, S., Hewson, C.M, Klaassen, R.H.G., Korner-Nievergelt, F., Kristensen, M.W & Hahn, S. (2012) Geolocation by light: Accuracy and precision affected by environmental factors. \emph{Methods in Ecology and Evolution}, doi: 10.1111/j.2041-210X.2012.00185.x

Wilson, R.P., Ducamp, J.J., Rees, G., Culik, B.M. & Niekamp, K. (1992) Estimation of location: global coverage using light intensity. \emph{Wildlife telemetry: remote monitoring and tracking of animals} (eds I.M. Priede & S.M. Swift), pp. 131-134. Ellis Horward, Chichester.
}
