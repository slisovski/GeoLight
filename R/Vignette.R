### Vignette

```{r}
library(GeoLight)
library(TwGeos)
```

### Download example dataset (source) and format to satisfy GeoLight format requirements

```{r}
# download.file('https://git.io/vrJgv', 'example_TAGS_format.csv')
exp <- read.csv("example_TAGS_format.csv", stringsAsFactors = F)
  exp$Twilight <- as.POSIXct(exp$datetime, tz = "UTC", 
                             format = "%Y-%m-%dT%T")
  exp$Rise <- ifelse(exp$twilight==1, TRUE, ifelse(exp$twilight==2, FALSE, NA))
exp <- subset(exp, !is.na(Rise), select = c("Twilight", "Rise"))  
  
  
twl <- export2GeoLight(exp)
```

### changeLight analysis

```{r}



```



  