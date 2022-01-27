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

myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT.Phillips = list(path_to_maxent.jar = getwd()))

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM'), #'MAXENT.Phillips'
  models.options = myBiomodOption,
  NbRunEval = 2,
  DataSplit = 70,
  Prevalence = 0.5,
  VarImport = 3,
  models.eval.meth = c('TSS','ROC'))

plot(myBiomodModelOut)

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
               
               
### Building ensemble-models
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.9))

# project models in space
myBiomodProjTime <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env = env_vars, 
                                      proj.name = "current",
                                      selected.models = 'all',
                                      binary.meth = 'TSS')

plot(myBiomodProjTime)

# ensemble-models projections on current time
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjTime )

plot(myBiomodEF)
