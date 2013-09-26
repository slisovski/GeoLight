\name{glfTrans}
\alias{glfTrans}
\title{Transformation of *.glf files}
\description{Transform *.glf files derived by the software GeoLocator for further analyses in \bold{\code{GeoLight}}.}
\usage{glfTrans(file)}

\arguments{
  \item{file}{the full patch and filename with suffix of the *.glf file.}
}
\value{A \code{data.frame} suitable for further use in \bold{\code{GeoLight}}.}

\details{The *.glf files produced by the software GeoLocator (Swiss Ornithological Institute) is a table with light intensity measurements over time. \code{glfTrans} produces a table with these measurements and transfer the data and time information into the format required by \bold{\code{GeoLight}} format (see: \code{\link{as.POSIXct}}).}
\seealso{\code{\link{gleTrans}}; \code{\link{luxTrans}} for transforming *.lux files produced by \emph{Migrate Technology Ltd}}
\author{
Simeon Lisovski
}