## Species distribution model ##
## Andre P. Silva ##
## January 26th, 2021 ##

# data setup
DataSpecies <- sdm_data #just changing name to match biomod script
myRespName <- "Prionailurus.bengalensis"
myResp <- as.numeric(DataSpecies[,myRespName])
head(myResp)
length(myResp)
unique(myResp)
myRespCoord = DataSpecies[,c("X_WGS84","Y_WGS84")]
head(myRespCoord)
nrow(myRespCoord)

# Initialisation
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = env_vars,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName)

#myBiomodOption <- BIOMOD_ModelingOptions(
#  MAXENT.Phillips = list(path_to_maxent.jar = getwd())) # to later add maxent.jar path

# Run individual models
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM'), #'MAXENT.Phillips'
  models.options = NULL, #myBiomodOption
  NbRunEval = 2,
  DataSplit = 70,
  Prevalence = 0.5,
  VarImport = 3,
  models.eval.meth = c('TSS','ROC'))

plot(myBiomodModelOut)

# Building ensemble-models
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.8),
  prob.mean = TRUE, # estimate the mean of probabilities
  prob.cv = TRUE, # coefficient variation
  prob.ci = TRUE, # confidence interval
  prob.ci.alpha = 0.05, # significance level for the confidence interval
  prob.median = TRUE, # estimate the median of probabilities
  prob.mean.weight = TRUE, # estimate the weighted mean of probabilities
  prob.mean.weight.decay = 'proportional') # bases weights on the model evaluation scores

# project models in current environment
myBiomodCurrent <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                     new.env = env_vars,
                                     proj.name = "current",
                                     selected.models = 'all',
                                     binary.meth = 'TSS')
#plot(myBiomodCurrent)

# ensemble-models projections on current environment
myBiomodEFCurrent <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                                projection.output = myBiomodCurrent)
plot(myBiomodEFCurrent)

# project models in future scenario Rcp26_70
myBiomodRcp26_70 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                     new.env = bio_Rcp26_70, 
                                     proj.name = "bio_Rcp26_70",
                                     selected.models = 'all',
                                     binary.meth = 'TSS')
plot(myBiomodRcp26_70)

myBiomodEFRcp26_70 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodRcp26_70)
plot(myBiomodEFRcp26_70)

# project models in future scenario Rcp26_70
myBiomodRcp85_70 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env = bio_Rcp85_70, 
                                      proj.name = "bio_Rcp85_70",
                                      selected.models = 'all',
                                      binary.meth = 'TSS')
plot(myBiomodRcp85_70)

myBiomodEFRcp85_70 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodRcp85_70)
plot(myBiomodEFRcp85_70)


# Model evaluation
## get AUC and TSS
get_evaluations(myBiomodModelOut)

## see what variables influence more AUC
get_variables_importance(myBiomodModelOut)

## species relationship with environmental variables
response.plot2(models = BIOMOD_LoadModels(myBiomodModelOut, models='GLM'),
               Data = get_formal_data(myBiomodModelOut,'expl.var'),
               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
               do.bivariate = FALSE,
               fixed.var.metric = 'median',
               col = brewer.pal(10, "Spectral"),
               legend = TRUE,
               data_species = get_formal_data(myBiomodModelOut,'resp.var'))



# plot ensemble forecasts
## load models manually 
list_projCurrent <- list.files("Prionailurus.bengalensis/proj_current/individual_projections", full.names = TRUE)
(list_proj <- grep("grd", list_proj, value=TRUE))
Pbcurrent <- raster::stack(list_proj)
names(Pbcurrent) <- grep("grd", list_proj, value=TRUE)
plot(Pbcurrent)

list_Rcp26_70  <- list.files("Prionailurus.bengalensis/proj_bio_Rcp26_70/individual_projections", full.names = TRUE)
(list_proj <- grep("grd", list_Rcp26_70, value=TRUE))
PbRcp26_70 <- raster::stack(list_proj)
names(PbRcp26_70) <- grep("grd", list_proj, value=TRUE)
plot(PbRcp26_70)

list_Rcp85_70  <- list.files("Prionailurus.bengalensis/proj_bio_Rcp85_70/individual_projections", full.names = TRUE)
(list_proj <- grep("grd", list_Rcp85_70, value=TRUE))
PbRcp85_70 <- raster::stack(list_proj)
names(PbRcp85_70) <- grep("grd", list_proj, value=TRUE)
plot(PbRcp85_70)

## join mean predictions and cv for the different future climatic scenarios in raster stack
Pbforecast <- raster::stack(Pbcurrent[[2]],
                            Pbcurrent[[1]],
                            PbRcp26_70[[2]],
                            PbRcp85_70[[2]])
names(Pbforecast) <- c("Pb_suitability_current",
                       "Pb_suitability_CV",
                       "Pb_suitability_Rcp26_70",
                       "Pb_suitability_Rcp85_70")
plot(Pbforecast) # typical final output map for your paper

