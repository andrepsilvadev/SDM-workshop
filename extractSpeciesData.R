## Species Data - builds input dataframe ##
## Filip Thorn and Andre P. Silva ##
## March 30th, 2020 ##

## extract GBIF data ------------------------------------------------------------

key <- name_backbone(name = "Mammalia")$usageKey # GBIF class key 
SEAsia<-c("CN;KH;BT;BD;IN;TH;MM;LA;VN;ID;MY;NP") #query in one dataframe
gbifData <- occ_search(taxonKey = key, # extract gbif mammal data
                       country = SEAsia,
                       hasCoordinate = TRUE,
                       year = '1970,2020',
                       limit=100000) # double check Records found by typing gbifData on console
#write_csv(gbifData, # used to build speciesDataFrame.R
#          file = "./data/gbifDataFull_20200918") # write as csv to avoid running this step all the time
#read gbifdataset (over 100MB > download from drive link, see README file in data folder)

## read datasets --------------------------------------------------------------

gbifData <- read_csv("./data/gbifDataFull_20200918.csv") # gbif
IUCN <- read.dbf("./data/TERRESTRIAL_MAMMALS.dbf", as.is = FALSE) # iucn this is the attribute table for the spatial polygon dataset avalible for terrestial mammals on IUCN
library(traitdata) # loads 22 datadets including pantheria and tetra_density
phylacine <- read_csv("./data/PHYLACINE_1.2.2/Traits/Trait_data.csv") # phylacine
edge <- read_csv("./data/EDScoresEDGE.csv") # evolutionary distinctiveness scores
zullinger <- read_csv("./data/Zullinger1984table3_edit.csv") # zullinger et al 1984
combined_imputed <- read_csv("./data/trait_data_imputed.csv") # COMBINE: a coalesced mammal database of intrinsic and extrinsic traits

## Prune GBIF data based on cutoff year and build input raster for rangeshifter with species initial
## locations within target landscape -removes duplicates within cells
spStack <- raster::stack()
data <- gbifData %>% dplyr::filter(year >= year) %>% drop_na(species)
n <- length(unique(data[["species"]]))
year <- 1970 # cut off year year to prune records
spStack <- raster::stack()
initDistrib <- data.frame(species = rep(NA,n),
                          family = rep(NA,n),
                          nlocations_noDuplic_afterCutoff = rep(NA,n)) 

initDistPath = "./data/spatial/initialDistribution"
for (i in 1:n){ 
  # Returns raster and dataframe with species records without duplicates per cell
  targetSpecies <- unique(data[["species"]])[[i]]
  sub <- data[grep(targetSpecies,data[["species"]]),]
  latlong <- data.frame(longitude = sub$decimalLongitude,
                        latitude = sub$decimalLatitude,
                        value = rep(1, nrow(sub)))
  spRaster <- rasterize(latlong[1:2], cellSize, field=latlong$value)
  names(spRaster) <- targetSpecies
  crop <- crop(spRaster, ext) # maybe it can be removed
  plot(crop, main = targetSpecies)
  rasterName <- paste(gsub(" ","_",targetSpecies), "InitDist.asc", sep="_")
  writeRaster(crop,
              filename=file.path(initDistPath, rasterName),
              format="ascii",
              NAflag = -9L) # requires changing cell size first cell value manually in txt...
  spStack <- addLayer(spStack, crop)
  initDistrib$species[i] <- as.character(targetSpecies)
  initDistrib$family[i] <- as.character(unique(sub$family))
  initDistrib$nlocations_noDuplic_afterCutoff[i] = sum(values(crop), na.rm=TRUE)
}

## filter species with ten or more records and plot
no.records <- 10
pruneLocat <- initDistrib %>% 
  dplyr::filter(nlocations_noDuplic_afterCutoff >= no.records) %>%
  arrange(desc(nlocations_noDuplic_afterCutoff))
(nrow(pruneLocat))
p <- ggplot(pruneLocat, aes(x=reorder(family,nlocations_noDuplic_afterCutoff), y=nlocations_noDuplic_afterCutoff)) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  scale_y_continuous(breaks = seq(0, 200, by = 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## join trait data from different sources --------------------------------------

spp <- pruneLocat %>% arrange(species)
spp <- spp$species

basedf <- tibble(
  index = 1:length(spp),
  species = spp,
  SpeciesCommon = rep(NA, length(species)))

iucnSub <- IUCN %>% rename(species = binomial) %>%
  distinct(species, .keep_all = TRUE) %>% 
  filter(species %in% spp) %>%
  arrange(species) %>%
  dplyr::select(species, order_, family, genus, category) %>% 
  mutate(order = tolower(order_),
         family = tolower(family),
         genus = tolower(genus), .keep = "unused")

pantheriaSub <- pantheria %>%
  unite("species", Genus:Species, sep = " ", remove = FALSE) %>%
  distinct(species, .keep_all = TRUE) %>% 
  filter(species %in% spp) %>%
  arrange(species)  %>%
  rename(order = Order,
         family = Family,
         genus = Genus)

tetra_densitySub <- tetra_density %>% # not used yet, to be used later for model validation
  filter(scientificNameStd %in% spp) %>%
  arrange(scientificNameStd)  %>%
  count(scientificNameStd)

phylacine$Binomial.1.2 <- gsub("_"," ",phylacine$Binomial.1.2)
phylacineSub <- phylacine %>% 
  rename(species = Binomial.1.2) %>%
  filter(species %in% spp) %>%
  arrange(species)

edgeSub <- edge %>% 
  filter(species %in% spp) %>%
  arrange(species)

zullingerSub <- zullinger %>% 
  rename(GompK = K,
         GompA = A,
         GompI = I) %>%
  filter(species %in% spp) %>%
  arrange(species)

joinedDatasets <- basedf %>% # dataframe with missing data
  left_join(iucnSub, by = "species", suffix = c("", ".iucn")) %>%
  left_join(pantheriaSub, by = "species", suffix = c("", ".pantheria")) %>%
  left_join(phylacineSub, by = "species", suffix = c("", ".phylacine")) %>%
  left_join(edgeSub, by = "species", suffix = c("", ".edge")) %>%
  left_join(zullingerSub, by = "species", suffix = c("", ".zullinger"))

## missing data inputation ------------------------------------------------------

sapply(joinedDatasets, function(x) sum(is.na(x))) # number NA per column
sapply(joinedDatasets, function(x) sum(is.na(x)/length(x))) # proportion NA per column

pantheriaMeans <- pantheriaSub %>%
  mutate(order = tolower(order),
         family = tolower(family),
         genus = tolower(genus)) %>% 
  group_by(genus) %>% # asumming phylogentic trait conservatism within genera
  summarise(meanMaxLong = mean(MaxLongevity_m, na.rm=T),
            nMaxLongGenera = sum(!is.na(MaxLongevity_m)),
            meanLitterSize = mean(LitterSize, na.rm=T),
            nMeanLitterSize = sum(!is.na(LitterSize)),
            meanLittersPerYear = mean(LittersPerYear, na.rm=T),
            nLittersPerYear = sum(!is.na(LittersPerYear)))

gompMeans <- zullinger %>% 
  mutate(order = tolower(order)) %>% 
  group_by(order) %>% # later this needs to be replaced by family
  summarise(meanGompK = mean(K),
            meanGompA = mean(A),
            meanGompI = mean(I),
            nGompGenera = sum(!is.na(K)))

inputNAdata <- joinedDatasets %>%
  mutate(AdultBodyMass_g = coalesce(AdultBodyMass_g, Mass.g)) %>% # complement Pantheria body mass with phylacine BM
  left_join(pantheriaMeans, by = "genus") %>% # fill Max Age NAs with genera means
  mutate(MaxLongevity_m = coalesce(MaxLongevity_m, meanMaxLong),
         LitterSize = coalesce(LitterSize, meanLitterSize),
         LittersPerYear = coalesce(LittersPerYear, meanLittersPerYear)) %>%
  left_join(gompMeans, by = "order") %>% # fill gomp NAs with order means
  mutate(GompK = coalesce(GompK, meanGompK),
         GompA = coalesce(GompA, meanGompA),
         GompI = coalesce(GompI, meanGompI))

sapply(test, function(x) sum(is.na(x))) # number NA per column
sapply(test, function(x) sum(is.na(x)/length(x))) # proportion NA per column

inputNAdata <- inputNAdata %>%
  drop_na(AdultBodyMass_g,
          MaxLongevity_m,
          meanGompK,
          LitterSize,
          LittersPerYear)

write_csv(inputNAdata, # used to build speciesDataFrame.R
          file = "./data/inputNAdata")
