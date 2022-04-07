future updates

### github minimal ----------------------------------------------------

create project
invite collaborators
download project to R
commit from R


### basic project folder structure ----------------------------------------------------

data
spatialData
scripts
  libraries.r
  environmentalData.r
  model.r
  pipeline.r
notebooks







### read multiple rasters and resample data -----------------------------------------------------

# create folder
folder <- "spatialData" # drop all the environmental variables in this folder
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}


readRaster <- function(folder, pattern){
  # function to read and stack all raster with the same pattern.
  # only rasters with same spatial extent and resolution can be stacked
  list.files <- list.files(path = folder,
                           pattern = pattern,
                           full.names = TRUE)
  stack <- raster::stack(list.files)
}


pattern <- "bio|elev"
bioclim <- readRaster(folder = folder,
                      pattern = pattern)

pattern <- "MODIS"
modis <- readRaster(folder = folder,
                    pattern = pattern)

modis <- projectRaster(modis, # change raster projection 
                       crs = crs(bioclim),
                       res = res(bioclim))

#extent <- raster::extent(-120, -30, -60, 30) # full extent #(xmin, xmax, ymin, ymax)
extent <- raster::extent(-70, -60, -10, -1) # test extent
bioclim <- raster::crop(bioclim, extent)
modis <- raster::crop(modis, extent)

# resample to bioclimatic raster resolution
modis <- resample(x = modis,
                  y = bioclim,
                  method = "ngb") # use ngb with integer data 

envVars <- raster::stack(bioclim, modis)
names(envVars) <- c("bio1","modis")
plot(envVars)


# check correlation between environmental variables
corr <- raster::layerStats(envVars, 'pearson', na.rm=T)
corr

# select non correlated variables
sel <- subset(env_vars, 1:2) # chose variable you want to keep