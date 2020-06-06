---
title: "Kriging with geoR"
author: "Agda Loureiro"
date: "5/13/2020"
output: html_document
runtime: shiny
---
# This code is made for geostatistical interpolation with ordinary kriging using geoR package

## For this we use the RStudio iteration with R version 3.6.3

#There are included the following analysis:

1st - First steps in R, installing and loading libraries, loading directory to source file location

2nd - Variogram analysis - we use Method of Moments (MoM) for experimental variogram with GLS fitting of theorical variogram with 
spheric, expanential and gaussian models.

3rd - Leave-one-out Cross validation (LOOCV) to choose semivariogram fit.

4th - Interpolation with kriging:  puntual Ordinary kriging (OK)

5th - Exporting interpolated maps

# 1. First steps in R:# 

## 1.1 - We start by cleaning R environment ##
```{r}
rm(list = ls())

gc(reset=T)

graphics.off()
```
## 1.2 - And install required packages
```{r}
#install.packages("pacman")

pacman::p_load(geoR, raster, rstudioapi, sp)
```

## 1.3 - Than we set working directory to source file location (directory turns to be the location where the R script is saved in your computer)

```{r}
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
```

## 1.4 - Loading data: our data is already free of outliers; we strongly recommend data preprocessing prior to interpolation

```{r}
data.n = read.csv(file = "../data/data points/data.csv", header = TRUE, sep = ',')
```

GeoR works wih geodata format, and it has a problem with high values of the spatial coordiantes. To solve this, we reduce the values of x and y by its minimum value
```{r}
data = data.frame(z = data.n$z,
                  x= data.n$x - min(data.n$x),
                  y = data.n$y - min(data.n$y))
```

## 1.5 - We separate the primary variable. This will facilitate analysis

```{r}
solo_atr<- "z"
```


## 1.6 - We transform data format to "geodata" 

```{r}
data.geo = as.geodata(data, coords.col = c(2,3),
                      data.col = solo_atr)
```

## 1.7 - Data visualization according to the "z" values

```{r}
plot(data.geo, lowess = TRUE)
points(data.geo)
```

# 2. Variogram analysis: using MoM.

## 2.1 - We start by visualization of the empirical variogram

```{r}
var_exp = variog(data.geo, max.dist = 1200)
plot(var_exp)
```
We can explore variogram by

vector containing distances:
```{r}
var_exp$u
```

vector with the semivariances values at u distances

```{r}
var_exp$v
```

number of pair of points by bins

```{r}
var_exp$n
```

### 2.1.1 If you have trend, use this code instead

For a first order trend, use trend argument as trend = "1st"
For a second order trend, use trend argument as trend = "2nd"

```{r}
var_exp = variog(data.geo, max.dist = 1200, trend = "1st")
plot(var_exp)
```


## 2.2 - Then we can eye-fit initial values for the semivariogram parameters. 

We make an eye-fit and save each value, it uses a an interactive Tcl-Tk interface. It will basecally open an interface to eye select values of semivariogram parameters.
 
```{r}
eye_fit_sph = eyefit(var_exp, silent = FALSE)
eye_fit_sph
```
 
```{r}
eye_fit_exp = eyefit(var_exp, silent = FALSE)
eye_fit_exp
```

```{r}
eye_fit_gau = eyefit(var_exp, silent = FALSE)
eye_fit_gau
```
sigmasq = partial sill
phi = range
tausq = nugget

## 2.3 - Experimental semivariogram by MoM and fitting of a parametric model to the semivariogram

We will test 4 models:

### 2.3.1 - Spheric

```{r}
fit_sph = geoR::variofit(var_exp, ini.cov.pars = eye_fit_sph, cov.model = "sph",
                         fix.nugget = F)
```

### 2.3.2 - Exponencial

```{r}
fit_exp = geoR::variofit(var_exp, ini.cov.pars = eye_fit_exp, fix.nug = F,
                         cov.model = "exp")
```

### 2.3.2 - Gaussian

```{r}
fit_gau = geoR::variofit(var_exp, ini.cov.pars = eye_fit_gau, cov.model = "gau",
                         fix.nugget = F)
```

# 3. LOOCV: leave-one-out cross validation to choose the best model.

### 3.1 - LOOCV of spheric model

```{r}
loocv_sph = xvalid(data.geo, model = fit_sph)
```

### 3.2 - LOOCV of exponential model

```{r}
loocv_exp = xvalid(data.geo, model = fit_exp)
```

### 3.3 - LOOCV of gaussian model

```{r}
loocv_gau = xvalid(data.geo, model = fit_gau)
```

### 3.4 Getting metrics

For Spheric

```{r}
metrics_sph = c(mae = hydroGOF::mae(loocv_sph$predicted, loocv_sph$data),
                         me = hydroGOF::me.data.frame(loocv_sph$predicted, loocv_sph$data),
                         rmse = hydroGOF::rmse(loocv_sph$predicted, loocv_sph$data),
r2 = hydroGOF::br2(loocv_sph$predicted, loocv_sph$data),
ave = hydroGOF::NSE(loocv_sph$predicted, loocv_sph$data,
                    ),willmott = hydroGOF::md(loocv_sph$predicted, loocv_sph$data))

```

For Exponential

```{r}
metrics_exp = c(mae = hydroGOF::mae(loocv_exp$predicted, loocv_exp$data),
                         me = hydroGOF::me.data.frame(loocv_exp$predicted, loocv_exp$data),
                         rmse = hydroGOF::rmse(loocv_exp$predicted, loocv_exp$data),
r2 = hydroGOF::br2(loocv_exp$predicted, loocv_exp$data),
ave = hydroGOF::NSE(loocv_exp$predicted, loocv_exp$data),
willmott = hydroGOF::md(loocv_exp
                        $predicted, loocv_exp$data))
```

For Gaussian

```{r}
metrics_gau = c(mae = hydroGOF::mae(loocv_gau$predicted, loocv_gau$data),
                         me = hydroGOF::me.data.frame(loocv_gau$predicted, loocv_gau$data),
                         rmse = hydroGOF::rmse(loocv_gau$predicted, loocv_gau$data),
r2 = hydroGOF::br2(loocv_gau$predicted, loocv_gau$data),
ave = hydroGOF::NSE(loocv_gau$predicted, loocv_gau$data),
willmott = hydroGOF::md(loocv_gau$predicted, loocv_gau$data))
```

And then adding all together

```{r}
metrics = data.frame(rbind(metrics_sph, metrics_exp, metrics_gau))
```

We can now export it

```{r}
write.csv(metrics, "../metrics/metrics.csv")
```

# 4 - Kriging interpolation

## 4.1 - We create a grid for interpolation

To performe this we open our data boundary/cotorno
```{r}
contorno <- shapefile("../data/contorno/cotorno.shp")
```

And then we create a grid

```{r}
r = raster::raster(contorno, res = 5) #  "res" sets pixel resolution

rp = raster::rasterize(contorno, r, 0) 

# convert raster to DataFrame
grid <- as.data.frame(x = rp, xy = T, na.rm = T)

summary(grid)  
```


change coords values for the same of the dt.geo

```{r}
grid$x <- grid$x - min(data$x)
grid$y <- grid$y - min(data$y) 
```

 collect the prediction grid into a geodata-object
 
```{r}
grid.geodata <- as.geodata(grid, data.col = 3,# ID
                           coords.col = c(1:2)) # coords (x and Y)

```
Ordinary Kriging (constant mean)
```{r}
ord.kriging <- krige.conv(geodata = data.geo, locations = grid.geodata$coords,
                          krige = krige.control(type.krige = "OK",
                                                obj.model = fit_exp)) # choose by resume analysis

kriging <- ord.kriging
```

Universal Kriging (trend, modelled by spatial coordinates)

```{r}
univ.kriging = geodata = krige.conv(data.geo, locations = grid.geodata$coords,
                          krige = krige.control(type.krige = "OK",
                                                obj.model = fit_exp,
                                                trend.d = ~ x + y,
                                                trend.l = ~ x + y)) 

```

input coord, predictions values into k.df

```{r}
k.df = data.frame(coords = grid[,1:2], predict = kriging$predict) 

colnames(k.df)[3] <- paste(solo_atr, sep = "_") 

```

 change coords values again

```{r}
k.df$coords.x <- k.df$coords.x + min(data$x)

k.df$coords.y <- k.df$coords.y + min(data$y)
```


# 5 - Export map


We first convert maps format to raster and add the maps projection

##5.1 - Convert to raster

```{r}
krig.map <- rasterFromXYZ (k.df) 
```


 set CRS

```{r}
crs(krig.map) <- crs(contorno)
```

 plot raster

```{r}
plot(krig.map, main = paste(solo_atr, "- kriging predictions", sep = " "))
```


## 5.1 - Exporting the map


```{r}
writeRaster(krig.map, "../mapas/Variavel_kriged",format="GTiff",
            overwrite = F, NAflag=-9999)
```


