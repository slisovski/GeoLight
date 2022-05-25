# Analysis of light based geolocator data

The GeoLight package provides basic functions for global positioning
based on light intensity measurements over time.  Positioning process
includes the determination of sun events, a discrimination of
residency and movement periods, the calibration of period-specific
data and, finally, the calculation of positions.

## Installing

The package is easily installed from GitHub, using the devtools package. 

```R
devtools::install_github("SLisovski/GeoLight")
```

If you don't have `devtools` installed already, install it first. 

```R
install.packages("devtools")
```

(GeoLight otherwise does not need devtools for normal use.)

## Known issues

The mergeSite2 functions seem to not be stable and may result in errors. The function is rather experimental and I do not recommend using it without a good look into the source code. 