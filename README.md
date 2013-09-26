# Analysis of light based geolocator data

The GeoLight package provides basic functions for global positioning
based on light intensity measurements over time.  Positioning process
includes the determination of sun events, a discrimination of
residency and movement periods, the calibration of period-specific
data and, finally, the calculation of positions.

## Installing

The package is easily built with RStudio

1. Install R

2. Install [RStudio](http://www.rstudio.com)

3. Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) or equivalent for your platform

4. Install [devtools](http://cran.r-project.org/web/packages/devtools/index.html) and [roxygen2](http://cran.r-project.org/web/packages/roxygen2/index.html) packages and dependencies in R

5. Install [maps](http://cran.r-project.org/web/packages/maps/index.html), [changepoint](http://cran.r-project.org/web/packages/changepoint/index.html) packages and dependencies.

6. Clone the repository from GitHub (https://github.com/mdsumner/GeoLight.git).

7. Create an Rstudio project in the folder containing this README file.

8. In the build tab, click `Configure Build Tools...` and click
`Generate documentation with Roxygen`, select `Configure` and choose to generate `Rd files` leaving the other options as they are.

9. Choose `Roxygenize` from the `Build` tab

10. Choose `Build & Reload` to make the package immediately available to R, or choose `More/Build source package` `More/Build binary package` from the `Build` tab to make source or binary packages.