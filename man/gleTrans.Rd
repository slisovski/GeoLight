\name{gleTrans}
\alias{gleTrans}
\title{Transformation of *.gle files}
\description{Function to transform *.gle files derived by the software GeoLocator for further analyses in \bold{\code{GeoLight}}.}
\usage{gleTrans(file)}

\arguments{
  \item{file}{the full patch and filename with suffix of the *.gle file.}
}
\value{A \code{data.frame} suitable for further use in \bold{\code{GeoLight}}.}

\details{The *.gle file derived by the software "GeoLocator" (Swiss Ornithological Institute) is a table with interpolated light intensities over time gathered from the *.glf file. Furthermore a column defines wether the light intensity passes the defined light intensity threshold in the morning (sunrise) or in the evening (sunset). This information is used in \code{gleTrans()}, to create a table with two subsequent twilight events (\emph{tFirst, tSecond}) and \emph{type} defining wether \emph{tFirst} refers to sunrise (1) or sunset (2). Date and time information will be transferred into a \code{GeoLight} appropriate format (see: \code{\link{as.POSIXct}}).}
\seealso{\code{\link{glfTrans}}}
\author{
Simeon Lisovski
}