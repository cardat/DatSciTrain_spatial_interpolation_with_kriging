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
library(rgdal)
library(sp)
library(automap)
library(raster)

######################################################################
infile_spatial <- "data_provided/nmmaps_city_subset_spatial_20220720.gpkg"
infile_source <- "data_provided/nmmaps_weather_city_subset_1987_2000_totals_20220720.csv"
outputfilename <- "data_derived/heatwave_1995_kriging_model.tif"
var <- "tmpd"
xname <- "long"
yname <- "lat"
city_i <- "Chicago"
year_i <- 1995
month_i <- 7
######################################################################

#### load the spatial data ####
## load the target study region
shp <- st_read(infile_spatial)

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

## visualise
with(input_data_csv[city == city_i & year == year_i & month == month_i,], plot(date, tmpd, type = "b"))

#### do: get avg  ####
indat <- input_data_csv[,.(tmpd = mean(tmpd, na.rm = T)), by = .(city, year, month)][year == year_i & month == month_i]

## prep the input spatial object to kriging
shp2 <- st_drop_geometry(shp)
shp2$city <- shp2$address
shp2$address <- NULL
indat2 <- merge(indat, shp2, by = "city")[rev(order(tmpd))]
indat2
## turn this from a data table into a data frame
setDF(indat2)
plot(indat2[,xname], indat2[,yname], col = heat.colors(nrow(indat2)), pch = 16)

# set up the projection info
epsg <- make_EPSG()
str(epsg)
# make spatial
input_data <- SpatialPointsDataFrame(coords = indat2[,c(xname, yname)],
                                     indat2, proj4string = CRS(epsg[epsg$code == 4326,"prj4"]))
str(input_data@data)
proj4string(input_data)
# NOTE there is something wrong with the proj4string.  ignore this today

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
bb <- st_bbox(shp)
bb
##plot(insp)
##axis(1); axis(2)
x.range <- c(bb[1] - 2, bb[3] + 2)
y.range <- c(bb[2] - 2, bb[4] + 2)
## CUSTOM EXTENT
##x.range <- c(150.25, 151.4) 
##y.range <- c(-34.3, -33.3)
## set resolution 0.25 degrees or about 25km
res <- 0.25
grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = res),
                   y = seq(from = y.range[1], to = y.range[2], by = res))  

coordinates(grd) <- ~x + y

gridded(grd) <- TRUE
str(grd)

proj4string(grd) <- CRS(proj4string(input_data))
### there is another warning about the projection prog4string problem.
## ignore this today

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
writeRaster(out_raster, outputfilename, format = "GTiff", overwrite=T)


#### Problems
## 1: create a figure in QGIS.
## 2: Add the values at the points using labels (hint: you can load a
## csv into QGIS if it has a lat and long)
