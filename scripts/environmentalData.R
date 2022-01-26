## Title: land data ##
## Author: Andre P. Silva ##
## Date: November 17th, 2020 ##

source("./scripts/R/customFunctions.R") # libraries

## extract annual land use rasters per variable -----------------------------------

var.luh <- c("primf", "primn", "secdf", "secdn", "pastr", "range",
             "urban", "c3ann", "c3per", "c4ann", "c4per", "c3nfx", "secma", "secmb") # variable definitions in annualScenarios/luh/LUH2_v2h_README.pdf
scenario.luh <- c("historical", "ssp245", "ssp585")
scenario.years <- list(850:2015, 2015:2100, 2015:2100)

## extract and save luh variable raster by year
for (i in 1:length(scenario.luh)) {
  for (j in 1:length(var.luh)) {
    year <- scenario.years[[i]]
    files <- list.files(luh.path)
    pattern <- scenario.luh[i]
    greprs <- grep(pattern, files, value=TRUE)
    rs <- raster::stack(file.path(luh.path, greprs), varname = var.luh[j])
    names(rs) <- paste(var.luh[j], year, sep="_")
    rs <- crop(rs, ext)
    if (pattern=="historical"){ # subsets luh historical raster to match climate data (1850 onward)
      rs <- rs[[1001:1166]]
    } else {
      rs <- rs
    }
    raster::writeRaster(rs,
                        filename= file.path("./data/spatial/luh/annual", paste0("luh_",scenario.luh [i], ".asc")),
                        bylayer=TRUE,
                        suffix=names(rs),
                        overwrite=TRUE)
  }
}

years <- 1850:2100
for (i in 1:length(years))  {
  grepYear <- grep(years[[i]], files, value = TRUE)
  luhYear <- stack(file.path("./data/spatial/luh/annual", grepYear))
  luhYear <- #stack(crop(luhYear, ext))
    luhYear.list[[i]] <- luhYear
  #SDM requires same name for calibration and projection variables - to check if layers have correct years comment names
  names(luhYear.list[[i]]) <- gsub("_[0-9]+","",names(luhYear.list[[i]]))
}


mergeLayerByName <- function(x, fun, indice, filename) {
  stackApply(x, indices, fun, filename, na.rm=TRUE)
}

EVsMerged <- lapply(EVs.projection.list, mergeLayerByName)


## climate scenarios --------------------------------------------------------------

# At the moment only done with one global climate model per scenario
# CMCC-CM was selected based on scores reported by McSweeney et al 2015 for Asia
#year <- c(1850, 1851, 1852) only to test
year <- seq(1850, 2100, by=1) # 1 year time steps
nyears <- length(year)
nMonths <- 12
nlayers <- 12*nyears 
indexYear <- seq(1, nlayers, by=nMonths)
scenario <- c("historical", "rcp45", "rcp85")
scenario <- scenario[1] # just for testing
var <- c("pr", "tasmax", "tasmin")
var <- var[1] # just for testing
#chelsa.path <- "./data/spatial/chelsa/"
#luh.path <- "./data/spatial/luh" 
#ext <- c(60,110,5,45) # SouthernAsiaExt - corresponds roughly to mammal pd regions 20,21,22,23,26 identified in Holt et al 2013
#ext <- c(96,100,26,29) # 192 cell test run
#cellSize <- raster(file.path(luh.path,"multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MESSAGE-ssp245-2-1-f_gn_2015-2100.nc"),
varname = "primf")[[1]] # resolution raster
#cellSize <- crop(cellSize, ext)
#ncell(cellSize)

##not sure if this correct NEEDs TO BE TESTED ! -> Juliano
##WARNING - before running doublecheck files order in source folder this will determine order they will be stacked
for (i in 1:length(var)) {
  path <- file.path(chelsa.path, var[i])
  for (j in 1:length(scenario)){ # read, stack, resample and extract annual climate data
    files <- list.files(path) # order of files is the order they will be stacked - files in folder need to be in correct order
    stack <- stackByPattern (scenario[j], files)
    stack <- crop(stack, ext)
    aggreg <- aggregate(stack, fact=5, fun=mean) # aggregate cells closer to luh resolution
    resamp <- resample(aggreg, cellSize, method="bilinear") #resample to exact luh resolution, bilinear continuous data
    mask <- mask(resamp, cellSize)
    annualMeanExtract(mask, nMonths, indexYear, paste0(scenario[j]))
  }
}

## quick test for clim loop using random raster - do not delete
#pr1852 <- stack(file.path(chelsa.path,"pr","CHELSAcmip5ts_pr_CMCC-CM_historical_1850-1869_V1.1.nc"))[[25:36]]
#pr1852 <- crop(pr1852, ext)
#aggreg1852 <- aggregate(pr1852, fact=5, fun=mean) # aggregate cells closer to luh resolution
#resamp1852 <- resample(aggreg1852, cellSize, method="bilinear")
#mask1852 <- mask(resamp1852, cellSize)
#mean_pr1852 <- mean(mask1852)
#plot(mean_pr1852)
#rasterToTest <- raster("./data/spatial/annualScenarios/pr_historical_annualMean_1852.asc") 
# Correlation between layers
#r1 <- mean_pr1852
#r2 <- rasterToTest
#cor(values(r1),
#    values(r2),
#    use = "na.or.complete")
#plot(stack(mask[[12]], masktest[[36]]))



## stack climate and luh rasters from target time steps ----------------------------

luhAnnual <- "./data/spatial/luh/annual"
files <- list.files(luhAnnual)

grepScenario <- grep("historical", files, value=TRUE)
targetYear <- as.vector(seq(1850, 2010, by=10))
historical <- lapply(targetYear, function(x) {
  grepYear <- grep(x, grepScenario, value=TRUE)
  stack <- raster::stack(file.path(luhAnnual, grepYear))
  return(stack)
})
names(historical) <- paste0("historical","_",targetYear)

grepScenario <- grep("ssp245", files, value=TRUE)
targetYear <- as.vector(seq(2020, 2090, by=10))
ssp245 <- lapply(targetYear, function(x) {
  grepYear <- grep(x, grepScenario, value=TRUE)
  stack <- raster::stack(file.path(luhAnnual, grepYear))
  return(stack)
})
names(ssp245) <- paste0("ssp245","_",targetYear)  

grepScenario <- grep("ssp585", files, value=TRUE)
targetYear <- as.vector(seq(2020, 2090, by=10))
ssp585 <- lapply(targetYear, function(x) {
  grepYear <- grep(x, grepScenario, value=TRUE)
  stack <- raster::stack(file.path(luhAnnual, grepYear))
  return(stack)
})
names(ssp585) <- paste0("ssp585", "_", targetYear) 

EVssp2 <- c(historical, ssp245)
EVssp5 <- c(historical, ssp585)







previous tests

#
prec <- raster::stack(prec.path)[[109]] #jan prec 2015
res <- resample(prec, primf[[1]], method="bilinear")
mask <- mask(res, primf[[1]])
plot(mask)

s <- stack(primf[[1]], secdf[[1]], mask) #year 2015
#names(s) <- c("primf","secdf")
plot(s)
(ncell(s)) #1036800
#Ext <- c(96,100,26,29) # 192 cell test run
Ext <- c(70,100,5,35) # SouthAsiaExt
EVs <- raster::stack(crop(s,Ext))

plot(EVs)
(ncell(EVs)) # 192; first trial with 6400 crashed
#LandscapeFile <- round(crp*100,2)
#writeRaster(round(crp,2), filename="./modelRuns/testSmallLandscape2/Inputs/landscape.asc",format="ascii",NAflag = -9L)

#open landscape file with wordpad and edit resolution to 25000 and first cell value to 0
#based on https://www.researchgate.net/post/How_to_convert_a_NetCDF4_file_to_GeoTIFF_using_R2
#https://luh.umd.edu/data.shtml 