#' @title do krige
#' @description a kriging model is useful to interpolate spatial data
#' from points to a surface. This is a Stretch Goal for the PUBH3007
#' Module 01 computer lab
#' @return a raster surface in R and saved to a geotiff file
#' @param infile_spatial file path to spatial data
#' @param infile_source file path to the input data
#' @param outputfilename file path to output geotiff file 
#' @param var variable name of column in mon_data 
#' @param xname  "long"
#' @param yname  "lat"
#' @param city_i choose one city to inspect visually
#' @param year_i choose the year to use for subset
#' @param month_i choose the month to use for subset

#### functions ####
library(data.table)
library(sf)
library(gstat)
library(sp)
library(automap)
library(raster)

######################################################################
infile_spatial <- "data_provided/air_pollution_monitoring_stations_by_city_20230917.csv"
infile_source <- "data_provided/air_pollution_monitoring_stations_dly_by_ste_australian_black_summer_20230917.csv"
outputfilename <- "data_derived/bushfire_smoke_australian_black_summer.tif"
var <- "pm25_final"
xname <- "lon"
yname <- "lat"
city_i <- "Sydney"
year_i <- 2019
month_i <- 12
day_i <- 10
date_i <- as.Date("2019-12-10")
######################################################################

#### load the spatial data ####
## load the target study region
shp <- fread(infile_spatial)
setDF(shp)
summary(shp)
shp$x <- shp$lon
shp$y <- shp$lat
shp <- st_as_sf(shp, coords = c("x", "y"), crs = 4326)

# load the attribute data
input_data_csv <- fread(infile_source)


#### clean and check ####
## inspect the data
plot(st_geometry(shp))
axis(1); axis(2)
shp
str(shp)
str(input_data_csv)

## create new variables
input_data_csv$year <- as.numeric(substr(input_data_csv$date, 1, 4))
input_data_csv$month <- as.numeric(substr(input_data_csv$date, 6, 7))


## prep the input spatial object to kriging
shp2 <- st_drop_geometry(shp)
shp2 <- shp2[shp2$city == "Sydney",]

indat <- input_data_csv[date == date_i]
indat2 <- merge(indat, shp2, by = c("station", "ste"))

## turn this from a data table into a data frame
setDF(indat2)
plot(indat2[,xname], indat2[,yname], col = heat.colors(nrow(indat2)), pch = 16)

# set up the projection info
# make spatial
input_data <- SpatialPointsDataFrame(coords = indat2[,c(xname, yname)],
                                     indat2, 
                                     proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs"))

str(input_data@data)
proj4string(input_data)

## Select data and make ready for input spatial
insp <- input_data[!is.na(input_data@data[,var]),]

#### do kriging ####
## This is how to autofit a variogram using the automap package
fmla <- reformulate('1', response = var)
vmod <- autofitVariogram(fmla, insp)
mod <- vmod$var_model

# Perform ordinary kriging, using the variogram calculated in “mod”
kmod <- gstat(NULL,var, fmla, insp, model=mod)

## For a raster prediction, load a template
## get the extent of the points
str(shp)
bb <- st_bbox(shp[shp$city == "Sydney",])
bb
##plot(insp)
##axis(1); axis(2)
x.range <- c(bb[1] - 0.2, bb[3] + 0.2)
y.range <- c(bb[2] - 0.2, bb[4] + 0.2)
## CUSTOM EXTENT
##x.range <- c(150.25, 151.4) 
##y.range <- c(-34.3, -33.3)
## set resolution 0.005 degrees or about 500m
res <- 0.005
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = res),
                   y = seq(from = y.range[1], to = y.range[2], by = res))  

coordinates(grd) <- ~x + y

gridded(grd) <- TRUE
str(grd)

proj4string(grd) <- CRS(proj4string(input_data))

## do a simple inverse distance weighted model and use this as the
## template grid for the kriging
idw <- idw(formula = fmla, locations = insp, newdata = grd, idp = 2)

r <- raster(idw)
out_raster <- interpolate(r, kmod)


#### do model checking ####
par(mfrow = c(1,2))   
plot(r)
plot(insp, add = T)
title("IDW")
plot(out_raster)
plot(insp, add = T)
title("Kriging")
## dev.off()


#### save the result ####
if(!dir.exists("data_derived")) dir.create("data_derived")
writeRaster(r, gsub(".tif", "_idw.tif", outputfilename), format = "GTiff", overwrite=T)
writeRaster(out_raster, gsub(".tif", "_krige.tif", outputfilename), format = "GTiff", overwrite=T)

#### Problems
## 1: create a figure in QGIS.
## 2: Add the values at the points using labels (hint: you can load a
## csv into QGIS if it has a lat and long)
