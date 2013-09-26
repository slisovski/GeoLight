i.sunelevation <-
function(lon, lat, year, month, day, hour, min, sec){

#-------------------------------------------------------------------------------
# lon: longitude in decimal coordinates
# lat: latitude in decimal coordinates
# year: numeric, e.g. 2006 (GMT)
# month: numeric, e.g. 8  (GMT)
# day: numeric, e.g. 6   (GMT)
# hour:  numeric e.g. 6  (GMT)
# min: numeric, e.g. 0   (GMT)
# sec: numeric, e.g.     (GMT)
#-------------------------------------------------------------------------------

datetime<-paste(year,"-", month,"-", day, " ", hour, ":", min, ":", sec, sep="")
gmt<-as.POSIXct(strptime(datetime, "%Y-%m-%d %H:%M:%S"), "UTC")
n <- gmt - as.POSIXct(strptime("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S"), "UTC")

# mean ecliptical length of sun 
L <- 280.46 + 0.9856474 * n
L <- as.numeric(L)

# Anomalie
g <- 357.528 + 0.9856003 * n
g <- as.numeric(g)

t.v <- floor(g/360)
g <- g - 360*t.v
g.rad <- g*pi/180

t.l <- floor(L/360)
L <- L - 360 * t.l
L.rad <- L*pi/180

# ecliptical length of sun
LAMBDA <- L + 1.915 * sin(g.rad) + 0.02*sin(2*g.rad)
LAMBDA.rad <- LAMBDA*pi/180

# coordinates of equator 
epsilon <- 23.439 - 0.0000004 * n
epsilon.rad <- as.numeric(epsilon)*pi/180

alpha.rad <- atan(cos(epsilon.rad)*sin(LAMBDA.rad)/cos(LAMBDA.rad))


alpha.rad <- ifelse(cos(LAMBDA.rad)<0, alpha.rad+pi, alpha.rad)
alpha <- alpha.rad*180/pi

deklination.rad <- asin(sin(epsilon.rad) * sin(LAMBDA.rad))
deklination <- deklination.rad*180/pi

# angle h
tag<-paste(year,"-", month,"-", day, " 00:00:00", sep="")
JD0<-as.POSIXct(strptime(tag, "%Y-%m-%d %H:%M:%S"), "GMT")
JD0 <- JD0 - as.POSIXct(strptime("2000-01-01 12:00:00", "%Y-%m-%d %H:%M:%S"), "GMT")
T0 <- JD0/36525

Time <- hour  + min/60  + sec/60/100
theta.Gh <- 6.697376 + 2400.05134 * T0 + 1.002738 * Time
theta.Gh <- as.numeric(theta.Gh)

t.d <- floor(theta.Gh/24)
theta.Gh <- theta.Gh-t.d*24

theta.G <- theta.Gh * 15

theta <- theta.G + lon      # Stundenwinkel des Frühlingspunktes
tau <- theta-alpha    # Stundenwinkel
tau.rad <- tau/180*pi

# Höhenwinkel h
h <- asin(cos(deklination.rad) * cos(tau.rad) * cos(lat/180*pi) + sin(deklination.rad) * sin(lat/180*pi))
h.grad <- h/pi*180

# correction because of refraction
R <- 1.02/(tan((h.grad+10.3/(h.grad+5.11))/180*pi))
hR.grad <- h.grad + R/60
return(hR.grad)
}

