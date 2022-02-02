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
myBiomodData_PArandom <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = env_vars,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName)

myBiomodData_PAdisk <-BIOMOD_FormatingData(resp.var,
                                    expl.var,
                                    resp.xy = NULL,
                                    resp.name = NULL,
                                    eval.resp.var = NULL,
                                    eval.expl.var = NULL,
                                    eval.resp.xy = NULL,
                                    PA.nb.rep = 0,
                                    PA.nb.absences = 1000,
                                    PA.strategy =c('random','disk'),
                                    PA.dist.min = 3000,
                                    PA.dist.max = 10000,
                                    PA.sre.quant = NULL,
                                    PA.table = NULL,
                                    na.rm = TRUE)

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
  eval.metric.quality.threshold = c(0.9),
  prob.mean = TRUE, # estimate the mean of probabilities
  prob.cv = TRUE, # coefficient variation
  prob.ci = TRUE, # confidence interval
  prob.ci.alpha = 0.05, # significance level for the confidence interval
  prob.median = TRUE, # estimate the median of probabilities
  prob.mean.weight = FALSE, # estimate the weighted mean of probabilities
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

BIOMOD_EnsembleForecasting( EM.output,
                            projection.output = NULL,
                            new.env = NULL,
                            xy.new.env = NULL,
                            selected.models ='all',
                            proj.name = NULL,
                            binary.meth = NULL,
                            filtered.meth = NULL,
                            compress = TRUE,...)

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




## Model evaluation
# get AUC and TSS
get_evaluations(myBiomodModelOut)

# see what variables influence more AUC
get_variables_importance(myBiomodModelOut)

# species relationship with environmental variables
response.plot2(models = BIOMOD_LoadModels(myBiomodModelOut, models='GLM'),
               Data = get_formal_data(myBiomodModelOut,'expl.var'),
               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
               do.bivariate = FALSE,
               fixed.var.metric = 'median',
               col = brewer.pal(10, "Spectral"),
               legend = TRUE,
               data_species = get_formal_data(myBiomodModelOut,'resp.var'))







## plot ensemble forecasts
# load models manually 
Pbcurrent <- raster("Prionailurus.bengalensis/proj_current/proj_current_Prionailurus.bengalensis_ensemble.grd")
PbRcp26_70 <- raster("Prionailurus.bengalensis/proj_bio_Rcp26_70/proj_bio_Rcp26_70_Prionailurus.bengalensis_ensemble.grd")
PbRcp85_70  <- raster("Prionailurus.bengalensis/proj_bio_Rcp85_70/proj_bio_Rcp85_70_Prionailurus.bengalensis_ensemble.grd")
# join predictions from the different scenarios in raster stack
Pbforecast <- raster::stack(Pbcurrent,PbRcp26_70,PbRcp85_70)
names(forecast) <- c("Pb_Ensemble_current","Pb_Ensemble_Rcp26_70","Pb_Ensemble_Rcp85_70")
plot(forecast)
