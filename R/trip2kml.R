#' Write a file which plots a trip in Google Earth
#' 
#' This function creates a .kml file from light intensity measurements over
#' time that can ve viewed as a trip in Google Earth.
#' 
#' 
#' @param file A character expression giving the whole path and the name of the
#' resulting output file including the .kml extension.
#' @param tFirst date and time of sunrise/sunset (e.g. 2008-12-01 08:30)
#' @param tSecond date and time of sunrise/sunset (e.g. 2008-12-01 17:30)
#' @param type either 1 or 2, defining \code{tFirst} as sunrise or sunset
#' respectively
#' @param degElevation sun elevation angle in degrees (e.g. -6 for "civil
#' twilight"). Either a single value, a \code{vector} with the same length as
#' \code{tFirst}.
#' @param col.scheme the color scheme used for the points. Possible color
#' schemes are: \code{\link{rainbow}}, \code{\link{heat.colors}},
#' \code{\link{topo.colors}}, \code{\link{terrain.colors}}.
#' @param point.alpha a \code{numerical value} indicating the transparency of
#' the point colors on a scale from 0 (transparent) to 1 (opaque).
#' @param cex \code{numerical value} for the size of the points.
#' @param line.col An character expression (any of \code{\link{colors}} or
#' hexadecimal notation), or numeric indicating the color of the line
#' connecting the point locations.
#' @return This function returns no data. It creates a .kml file in the in the
#' defined path.
#' @author Simeon Lisovski and Michael U. Kemp
#' @examples
#' 
#' data(hoopoe2)
#' attach(hoopoe2)
#' filter <- distanceFilter(tFirst,tSecond,type,distance=30)
#' trip2kml("trip.kml", tFirst[filter], tSecond[filter], type[filter], 
#' 		degElevation=-6, col.scheme="heat.colors", cex=0.7, 
#' 		line.col="goldenrod")
#' 
#' @export trip2kml
trip2kml <- function(file, tFirst, tSecond, type, degElevation, col.scheme="heat.colors", point.alpha=0.7, cex=1, line.col="goldenrod")
{
	if((length(tFirst)+length(type))!=(length(tSecond)+length(type))) stop("tFirst, tSecond and type must have the same length.")
	
	coord   <- coord(tFirst,tSecond,type,degElevation,note=F)
		index <- !is.na(coord[,2])
	datetime <- as.POSIXct(strptime(paste(ifelse(type==1,substring(tFirst,1,10),substring(tSecond,1,10)),
				" ",ifelse(type==1,"12:00:00","00:00:00"),sep=""),format="%Y-%m-%d %H:%M:%S"),"UTC")
		
	coord   <- coord[index,]
	longitude <- coord[,1]
	latitude <- coord[,2]

	date <- unlist(strsplit(as.character(datetime[index]), split = " "))[seq(1, 
        ((length(datetime[index]) * 2) - 1), by = 2)]
    time <- unlist(strsplit(as.character(datetime[index]), split = " "))[seq(2, 
        ((length(datetime[index]) * 2)), by = 2)]
	
	if(length(!is.na(coord[,2]))<1) stop("Calculation of coordinates results in zero spatial information.")
	
	if ((col.scheme%in% c("rainbow", "heat.colors", "terrain.colors", "topo.colors", 
						  "cm.colors"))==F) stop("The col.scheme has been misspecified.")
	
	seq   <- seq(as.POSIXct(datetime[1]),as.POSIXct(datetime[length(datetime)]),by=12*60*60)
	index2<- ifelse(!is.na(merge(data.frame(d=datetime[index],t=TRUE),data.frame(d=seq,t=FALSE),by="d",all.y=T)[,2]),TRUE,FALSE)
	
	usable.colors <- strsplit(eval(parse(text = paste(col.scheme, 
            "(ifelse(length(index2) < 1000, length(index2), 1000), alpha=point.alpha)", 
            sep = ""))), split = "")[index2]
	
	usable.line.color <- strsplit(rgb(col2rgb(line.col)[1, 
        1], col2rgb(line.col)[2, 1], col2rgb(line.col)[3, 
        1], col2rgb(line.col, alpha = 1)[4, 1], maxColorValue = 255), 
        split = "")
	
	date <- unlist(strsplit(as.character(datetime), split = " "))[seq(1, 
        ((length(datetime) * 2) - 1), by = 2)]
    time <- unlist(strsplit(as.character(datetime), split = " "))[seq(2, 
        ((length(datetime) * 2)), by = 2)]
    scaling.parameter <- rep(cex, length(latitude))
    
    data.variables <- NULL
    filename <- file
    write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>", filename)
    write("<kml xmlns=\"http://www.opengis.net/kml/2.2\">", filename, 
        append = TRUE)
    write("<Document>", filename, append = TRUE)
    write(paste("<name>", filename, "</name>", sep = " "), 
        filename, append = TRUE)
    
    write("  <open>1</open>", filename, append = TRUE)
    write("\t<description>", filename, append = TRUE)
    write("\t  <![CDATA[Generated using <a href=\"http://simeonlisovski.wordpress.com/geolight\">GeoLight</a>]]>", 
        filename, append = TRUE)
    write("\t</description>", filename, append = TRUE)
    write("<Folder>", filename, append = TRUE)
    write("  <name>Points</name>", filename, append = TRUE)
    write("<open>0</open>", filename, append = TRUE)
    for (i in 1:length(latitude)) {
        write("<Placemark id='point'>", filename, append = TRUE)
        write(paste("<name>", as.character(as.Date(datetime[i])), "</name>", sep = ""), 
            filename, append = TRUE)
        write("  <TimeSpan>", filename, append = TRUE)
        write(paste("    <begin>", date[i], "T", time[i], "Z</begin>", 
            sep = ""), filename, append = TRUE)
        write(paste("    <end>", date[ifelse(i == length(latitude), 
            i, i + 1)], "T", time[ifelse(i == length(latitude), 
            i, i + 1)], "Z</end>", sep = ""), filename, append = TRUE)
        write("  </TimeSpan>", filename, append = TRUE)
        write("<visibility>1</visibility>", filename, append = TRUE)
        write("<description>", filename, append = TRUE)
        write(paste("<![CDATA[<TABLE border='1'><TR><TD><B>Variable</B></TD><TD><B>Value</B></TD></TR><TR><TD>Date/Time</TD><TD>", 
                datetime[i], "</TD></TR><TR><TD>lat long</TD><TD>", 
                paste(latitude[i], longitude[i], sep = " "), 
                "</TABLE>]]>", sep = "", collapse = ""), filename, 
                append = TRUE)
        write("</description>", filename, append = TRUE)
        write("\t<Style>", filename, append = TRUE)
        write("\t<IconStyle>", filename, append = TRUE)
        write(paste("\t\t<color>", paste(noquote(usable.colors[[i]][c(8, 
            9, 6, 7, 4, 5, 2, 3)]), collapse = ""), "</color>", 
            sep = ""), filename, append = TRUE)
        write(paste("  <scale>", scaling.parameter[i], "</scale>", 
            sep = ""), filename, append = TRUE)
        write("\t<Icon>", filename, append = TRUE)
        write("\t\t<href>http://maps.google.com/mapfiles/kml/pal2/icon26.png</href>", 
            filename, append = TRUE)
        write("\t</Icon>", filename, append = TRUE)
        write("\t</IconStyle>", filename, append = TRUE)
        write("\t</Style>", filename, append = TRUE)
        write("\t<Point>", filename, append = TRUE)
        write(paste("\t<altitudeMode>", "relativeToGround", "</altitudeMode>", 
            sep = ""), filename, append = TRUE)
        write("<tesselate>1</tesselate>", filename, append = TRUE)
        write("<extrude>1</extrude>", filename, append = TRUE)
        write(paste("\t  <coordinates>", longitude[i], ",", latitude[i], 
            ",", 1, 
            "</coordinates>", sep = ""), filename, append = TRUE)
        write("\t</Point>", filename, append = TRUE)
        write(" </Placemark>", filename, append = TRUE)
    }
   
    write("</Folder>", filename, append = TRUE)
    write("<Placemark>", filename, append = TRUE)
    write("  <name>Line Path</name>", filename, append = TRUE)
    write("  <Style>", filename, append = TRUE)
    write("    <LineStyle>", filename, append = TRUE)
    write(paste("\t<color>", paste(noquote(usable.line.color[[1]][c(8, 
        9, 6, 7, 4, 5, 2, 3)]), collapse = ""), "</color>", sep = ""), 
        filename, append = TRUE)
    
    
    write(paste("      <width>1</width>", sep = ""), filename, 
        append = TRUE)
    write("    </LineStyle>", filename, append = TRUE)
    write("  </Style>", filename, append = TRUE)
    write("  <LineString>", filename, append = TRUE)
    write("    <extrude>0</extrude>", filename, append = TRUE)
    write("    <tessellate>1</tessellate>", filename, append = TRUE)
    write(paste("\t<altitudeMode>clampToGround</altitudeMode>", 
        sep = ""), filename, append = TRUE)
    write(paste("     <coordinates>", noquote(paste(longitude, 
        ",", latitude, sep = "", collapse = " ")), "</coordinates>", 
        sep = ""), filename, append = TRUE)
    write("    </LineString>", filename, append = TRUE)
    write("</Placemark>", filename, append = TRUE)
    write("</Document>", filename, append = TRUE)
    write("</kml>", filename, append = TRUE)
}
