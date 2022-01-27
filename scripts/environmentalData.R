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
env_vars <- raster::stack(env_vars)
env_vars

# join species and environmental data
targetspecies <- unique(gbifData$data$species)
emptystack <- raster::stack()
for(i in 1:length(targetspecies)){
  subset <- gbifData$data %>% dplyr::filter(species == targetspecies[i])
  xy <- data.frame(x=subset$decimalLongitude,
                   y=subset$decimalLatitude)
  spRaster <- rasterize(xy, env_vars, fun='count')
  emptystack <- addLayer(emptystack, spRaster)
}
plot(emptystack)
names(emptystack) <- unique(gbifData$data$species)
plot(emptystack)
# reclassify species raster
m <- c(NA, NA, 0, 0, +Inf, 1)
(rclmat <- matrix(m, ncol=3, byrow=TRUE)) # criteria for reclassification
rc <- reclassify(emptystack, rclmat)
plot(rc)
spStack <- rc # just changing stack name

full_data <- stack(spStack, env_vars)
plot(full_data)

xylandscape <- raster::coordinates(env_vars) 
colnames(xylandscape) <- c("X_WGS84","Y_WGS84") # to math with biomod script
sdm_data <- raster::extract(full_data, # raster or rasterstack
                            xylandscape, # landscape coordinates
                            method='simple', # or "bilinear" - value of the four nearest raster cells
                            buffer=NULL, # in meters (long lat) or map_units
                            small=FALSE,
                            cellnumbers=TRUE,
                            na.rm=TRUE,
                            df=TRUE) # return dataframe
sdm_data
sdm_data <- cbind(xylandscape, sdm_data) #join cell coordinates
head(sdm_data)
nrow(sdm_data)
sdm_data <- sdm_data %>% tidyr::drop_na()
sdm_data
nrow(sdm_data)
# not needed juts to double check
unique(sdm_data$Prionailurus.bengalensis)
unique(sdm_data$Prionailurus.rubiginosus)
unique(sdm_data$Prionailurus.viverrinus)
unique(sdm_data$Prionailurus.planiceps)
sum(sdm_data$Prionailurus.bengalensis)
sum(sdm_data$Prionailurus.rubiginosus)
sum(sdm_data$Prionailurus.viverrinus)
sum(sdm_data$Prionailurus.planiceps)

## exploring a bit to get to know the data set
#ggplot(full_data, aes(x=bio10, y=cells, colour=species)) +
#  geom_point()


