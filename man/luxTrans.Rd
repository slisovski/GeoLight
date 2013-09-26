\name{luxTrans}
\alias{luxTrans}
\title{Transformation of *.lux files}
\description{Transform *.lux files derived from \emph{Migrate Technology Ltd} geolocator deviced for further analyses in \bold{\code{GeoLight}}.}
\usage{luxTrans(file)}

\arguments{
  \item{file}{the full patch and filename with suffix of the *.lux file.}
}
\value{A \code{data.frame} suitable for further use in \bold{\code{GeoLight}}.}

\details{The *.lux files produced by \emph{Migrate Technology Ltd} are table with light intensity measurements over time. \code{luxTrans} produces a table with these measurements and transfer the data and time information into the format required by \bold{\code{GeoLight}} format (see: \code{\link{as.POSIXct}}).}
\seealso{\code{\link{gleTrans}} for transforming *.glf files produced by the software GeoLocator (\emph{Swiss Ornithological Institute})}
\author{
Simeon Lisovski
}