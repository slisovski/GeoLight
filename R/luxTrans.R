#' Transformation of *.lux files
#' 
#' Transform *.lux files derived from \emph{Migrate Technology Ltd} geolocator
#' deviced for further analyses in \bold{\code{GeoLight}}.
#' 
#' The *.lux files produced by \emph{Migrate Technology Ltd} are table with
#' light intensity measurements over time. \code{luxTrans} produces a table
#' with these measurements and transfer the data and time information into the
#' format required by \bold{\code{GeoLight}} format (see:
#' \code{\link{as.POSIXct}}).
#' 
#' @param file the full patch and filename with suffix of the *.lux file.
#' @return A \code{data.frame} suitable for further use in
#' \bold{\code{GeoLight}}.
#' @author Simeon Lisovski
#' @seealso \code{\link{gleTrans}} for transforming *.glf files produced by the
#' software GeoLocator (\emph{Swiss Ornithological Institute})
#' @export luxTrans
luxTrans <-
function(file="/path/file.lux") {
	
lux1 <- read.table(file,sep="\t",skip=21,col.names=c("datetime","time")) # read file
lux <- data.frame(datetime=as.POSIXct(strptime(lux1[,1],format="%d/%m/%Y %H:%M:%S",tz="UTC")),light=lux1[,2])

return(lux)
}

