# title do krige
# description a kriging model is useful to interpolate spatial data
#   from points to a surface. 
# returns a raster surface in R and saved to a geotiff file

# other notes:
# infile_spatial file path to spatial data
# infile_source file path to the input data
# outputfilename file path to output geotiff file 
# var variable name of column in mon_data 
# xname  "long"
# yname  "lat"
# city_i choose one city to inspect visually
# year_i choose the year to use for subset
# month_i choose the month to use for subset


#### functions ####
library(data.table)
library(terra)
library(sf)
library(stars)
library(gstat)
library(automap)

######################################################################
infile <- "data_provided/napmd_pm25_hrly_avg_NSW_20191210_0200_black_summer.csv"
outputfilename <- "data_derived/bushfire_smoke_australian_black_summer.tif"

value <- "value"
xname <- "lon"
yname <- "lat"

# define variable of interest
response_var <- "pm25"

# define study extent
# Sydney region
study_bounds <- c(xmin = 149.4719, ymin = -34.83117, 
                  xmax = 152.1305, ymax = -32.49607)
year_i <- 2019
month_i <- 12
day_i <- 10
hour_i <- 2
datetime_i <- as.POSIXct(paste0(year_i, "-", month_i, "-", day_i, " ", hour_i, ":00"),
                         tz = "UTC")

######################################################################

#### load the data ####
dat <- fread(infile)
sf.mons <- st_as_sf(dat, coords = c(xname, yname), crs = 4283)

#### clean and check ####
## inspect the data
plot(sf.mons["station"],
     pch = 20, col = "black",
     axes = T)
sf.mons

str(sf.mons)


## Select data and make ready for input spatial
insp <- sf.mons[(sf.mons$variable == response_var &
                   sf.mons$date_time_utc == datetime_i &
                   !is.na(value)), ]
insp <- st_crop(insp, st_bbox(study_bounds, crs = 4283))
plot(insp["value"], pch = 16, axes = T)



#### do kriging ####
## This is how to autofit a variogram using the automap package, requires use of SpatialPointsDataFrame
insp <- as(insp, "Spatial")
fmla <- reformulate('1', response = value)
vmod <- autofitVariogram(formula = fmla, 
                         input_data = insp)
mod <- vmod$var_model

# Perform ordinary kriging, using the variogram calculated in “mod”
kmod <- gstat(NULL, id = value, fmla, insp, model=mod)



#### Predict on grid ####
# create grid template
res <- 0.005
grd <- stars::st_as_stars(st_bbox(study_bounds), dx = res, dy = res)
grd <- st_set_crs(grd, 4283)

## do a simple inverse distance weighted model
r_idw <- idw(formula = fmla, locations = insp, newdata = grd, idp = 2)

## do the kriging model raster
r_krige <- predict(kmod, grd)



#### do model checking ####
# plot with terra
r_idw.terra <- rast(r_idw["var1.pred"])
r_krige.terra <- rast(r_krige["value.pred"])
v_insp.terra <- vect(insp)

par(mfrow = c(1, 2))
plot(r_idw.terra, col = rev(map.pal("ryg")), main = "IDW")
plot(v_insp.terra, "value", pch = 20, legend = F, col = "black", add = T)

plot(r_krige.terra, col = rev(map.pal("ryg")), main = "Kriging")
plot(v_insp.terra, "value", pch = 20, legend = F, col = "black", add = T)



#### save the result ####
if(!dir.exists("data_derived")) dir.create("data_derived")
writeRaster(r_idw.terra, gsub(".tif", "_idw.tif", outputfilename), overwrite=T)
writeRaster(r_krige.terra, gsub(".tif", "_krige.tif", outputfilename), overwrite=T)



#### Problems
## 1: create a figure in QGIS.
## 2: Add the values at the points using labels (hint: you can load a
## csv into QGIS if it has a lat and long)
