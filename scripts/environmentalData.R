## Title: environmentalData ##
## Author: Andre P. Silva ##
## Date: January 26th, 2022 ##


r <- raster::getData("worldclim", var="bio", res=10) # what is the most adequate resolution? #resolution is in degrees
plot(r)
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
clim_vars <- c("bio10","bio11","bio12") # How many variables shall I use?
subset_clim <- r[[clim_vars]] # subset target climatic variables 
plot(subset_clim)

#reduce extent to target region. But how to define the target region?
extent <- raster::extent(50, 140, -10, 40) #(xmin, xmax, ymin, ymax)
r_crop <- raster::crop(subset_clim, extent)
plot(r_crop)

b1 <- Raster::stackApply(r_crop, indices=c(1,2), fun=function(x){x/10}) # error here
b1 <- raster::stackApply(s, indices=c(1,1,1,2,2,2), fun=sum)

# check correlation between environmental variables
corr <- raster::layerStats(r_crop, 'pearson', na.rm=T)
env_vars <- raster::dropLayer(r_crop, "bio11")
env_vars


# this step is not necessary to run the SDMs - other analyses or exploratory
xy <- data.frame(x=gbifData$data$decimalLongitude,
                 y=gbifData$data$decimalLatitude)
head(xy)

extract_data <- raster::extract(env_vars, # raster or rasterstack
                            xy, # coordinates
                            method='simple', # or "bilinear" - value of the four nearest raster cells
                            buffer=NULL, # in meters (long lat) or map_units
                            small=FALSE,
                            cellnumbers=TRUE,
                            na.rm=TRUE,
                            df=TRUE) # return dataframe
head(extract_data)
nrow(extract_data)

full_data <- data.frame(species = gbifData$data$species,
                       lat = gbifData$data$decimalLatitude,
                       long = gbifData$data$decimalLongitude,
                       cells = extract_data$cells ,
                       bio10 = extract_data$bio10,
                       bio12 = extract_data$bio12)

head(full_data)
# exploring a bit to get to know the data set
ggplot(full_data, aes(x=bio10, y=cells, colour=species)) +
  geom_point()

