## Extract species data ##
## Andre P. Silva ##
## January 26th, 2022 ##

## extract GBIF data ------------------------------------------------------------
(key <- rgbif::name_backbone(name = "Prionailurus")$usageKey) # GBIF class key 

# data filtering
region <- c("CN;KH;BT;BD;IN;TH;MM;LA;VN;ID;MY;NP")
date <- '1970,2020'
# search data
gbifData <- rgbif::occ_search(taxonKey = key,
                       country = region,
                       hasCoordinate = TRUE,
                       year = date,
                       limit = 100000)
gbifData$data
#readr::write_csv(gbifData$data,
#          file = "./data/gbifData_20220126") # write as csv to avoid running this step all the time
#gbifData <- readr::read_csv("./data/gbifDataFull_20200918.csv") # gbif

# alternatively you may want to load the data directly from your files
#cats <- read_csv("./data/cats_RD.csv")
#head(cats)###################### ERROR HERE 

# visualize number of locations
(p <- ggplot(gbifData$data, aes(x=species)) +
  geom_histogram(stat="count"))

# plot in world map
worldmap <- rworldmap::getMap(resolution = "coarse")
plot(worldmap, col = "lightgrey", 
     fill = T, border = "darkgray",
     xlim = c(100, 140), ylim = c(-10,40), #xlim = c(80, 180), ylim = c(-50, 90),
     bg = "aliceblue")
maps::map.axes()
points(gbifData$data$decimalLongitude, gbifData$data$decimalLatitude, 
       col = brewer.pal(4,"Spectral"), cex = 1)
#legend(120, 30, legend=unique(gbifData$data$species)) #not working


