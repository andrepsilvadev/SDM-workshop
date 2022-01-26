## Species distribution model ##
## Andre P. Silva ##
## January 26th, 2021 ##

needs another format landscape zeros AND ONES

DataSpecies <- full_data[c("species", "lat", "long")]
head(DataSpecies)

myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[,c("X_WGS84","Y_WGS84")]
### Initialisation
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = EVs.calibration,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName)

myBiomodOption <- BIOMOD_ModelingOptions(
  #MAXENT.Phillips = list(path_to_maxent.jar = getwd()))
  MAXENT.Phillips = list(path_to_maxent.jar = 'Z:/MechSpatCons/data/spatial/SDMoutput/maxent.jar'))

myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'RF','GBM'), # "MAXENT.Phillips"
  models.options = myBiomodOption,
  NbRunEval = 10,
  DataSplit = 70,
  Prevalence = 0.5,
  VarImport = 3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE, # see what this is about
  do.full.models = FALSE,
  modeling.id = paste(myRespName, "FirstModeling", sep=""))

get_evaluations(myBiomodModelOut)

get_variables_importance(myBiomodModelOut)

response.plot(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE)

### Building ensemble-models
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7), 
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional')

# Make ensemble-models projections on current variable
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjTime,
  new.env = NULL,
  xy.new.env = NULL,
  selected.models = 'all',
  proj.name = NULL,
  binary.meth = NULL,
  filtered.meth = NULL,
  compress = TRUE)








MyBiomodSF <- function(sp.names){
  #setwd(SDMpath)
  myRespName = gsub(" ",".", sp.names)
  cat('\n', myRespName, 'modeling...')
  
  ### definition of data
  ## i.e keep only the column of our species
  myResp <- as.numeric(DataSpecies[,myRespName])
  myRespCoord = DataSpecies[,c("X_WGS84","Y_WGS84")]
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = EVs.calibration,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName)
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions(
    #MAXENT.Phillips = list(path_to_maxent.jar = getwd()))
    MAXENT.Phillips = list(path_to_maxent.jar = 'Z:/MechSpatCons/data/spatial/SDMoutput/maxent.jar'))
  #file.path(SDMpath,"maxent.jar")
  ### Modelling
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('GLM', 'RF','GBM'), # "MAXENT.Phillips"
    models.options = myBiomodOption,
    NbRunEval = 10,
    DataSplit = 70,
    Prevalence = 0.5,
    VarImport = 3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = FALSE, # see what this is about
    do.full.models = FALSE,
    modeling.id = paste(myRespName, "FirstModeling", sep=""))
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
                 file=file.path(myRespName,
                                paste(myRespName,
                                      "_formal_models_evaluation.txt", sep="")))
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(myRespName,
                                paste(myRespName,
                                      "_formal_models_variables_importance.txt", sep="")))
  
  ### Building ensemble-models
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = c(0.6), #tutorials often use 0.7
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional')
  
  ### Projection on environmental conditions over time
  for (i in 1:length(EVs.projection.list))  {
    myBiomodProjTime <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = EVs.projection.list[[i]], #it needs the same variable names as the calibration variables!
                                          proj.name = names(EVs.projection.list)[[i]],
                                          selected.models = 'all',
                                          binary.meth = 'TSS',
                                          compress = FALSE,
                                          build.clamping.mask = TRUE)
    
    # Make ensemble-models projections on current variable
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProjTime,
      new.env = NULL,
      xy.new.env = NULL,
      selected.models = 'all',
      proj.name = NULL,
      binary.meth = NULL,
      filtered.meth = NULL,
      compress = TRUE)
  }
}

sp.names <- spData$Species
system.time(biomod_multispecies <- lapply(sp.names, MyBiomodSF))














year.calibration <- 2000
EVs.calibration <- historical[["historical_2000"]]
var.names <- c("c3ann", "c3nfx", "c3per", "c4ann", "c4per", "pastr", "primf", # SDM requires same name for calibration and projection variables - to check if layers have correct years comment names
               "primn", "range", "secdf", "secdn", "secma", "secmb", "urban") #needs to follow this exact order to match order from grep - way to dit automatic potential problems
names(EVs.calibration) <- var.names

EVs.projection.list <- c(historical, ssp245) # needs to add ssp585
EVs.projection.list <- lapply(EVs.projection.list, function(x) {
  names(x) <- var.names
  return(x)})
spStack <- stack(list.files("./data/spatial/initialDistribution", full.names = TRUE))

coordinates <- coordinates(EVs.calibration[[1]])
rasValue <- raster::extract(spStack, coordinates)
DataSpecies <- cbind(coordinates, rasValue)
DataSpecies[is.na(DataSpecies)] <- 0
colnames(DataSpecies) <- gsub("_InitDist", "", colnames(DataSpecies))
colnames(DataSpecies) <-  gsub("_", ".", colnames(DataSpecies))
colnames(DataSpecies)[1:2] <- c('X_WGS84','Y_WGS84')

MyBiomodSF <- function(sp.names){
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(SDMpath)
  #myRespName = sp.names
  myRespName = gsub(" ",".", sp.names)
  cat('\n', myRespName, 'modeling...')
  
  ### definition of data
  ## i.e keep only the column of our species
  myResp <- as.numeric(DataSpecies[,myRespName])
  myRespCoord = DataSpecies[,c("X_WGS84","Y_WGS84")]
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = EVs.calibration,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName)
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions(
    #MAXENT.Phillips = list(path_to_maxent.jar = getwd()))
    MAXENT.Phillips = list(path_to_maxent.jar = 'Z:/MechSpatCons/data/spatial/SDMoutput/maxent.jar'))
  #file.path(SDMpath,"maxent.jar")
  ### Modelling
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('GLM', 'RF','GBM'), # "MAXENT.Phillips"
    models.options = myBiomodOption,
    NbRunEval = 10,
    DataSplit = 70,
    Prevalence = 0.5,
    VarImport = 3,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = FALSE, # see what this is about
    do.full.models = FALSE,
    modeling.id = paste(myRespName, "FirstModeling", sep=""))
  ### save models evaluation scores and variables importance on hard drive
  capture.output(get_evaluations(myBiomodModelOut),
                 file=file.path(myRespName,
                                paste(myRespName,
                                      "_formal_models_evaluation.txt", sep="")))
  capture.output(get_variables_importance(myBiomodModelOut),
                 file=file.path(myRespName,
                                paste(myRespName,
                                      "_formal_models_variables_importance.txt", sep="")))
  
  ### Building ensemble-models
  myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('TSS'),
    eval.metric.quality.threshold = c(0.6), #tutorials often use 0.7
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional')
  
  ### Projection on environmental conditions over time
  for (i in 1:length(EVs.projection.list))  {
    myBiomodProjTime <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = EVs.projection.list[[i]], #it needs the same variable names as the calibration variables!
                                          proj.name = names(EVs.projection.list)[[i]],
                                          selected.models = 'all',
                                          binary.meth = 'TSS',
                                          compress = FALSE,
                                          build.clamping.mask = TRUE)
    
    # Make ensemble-models projections on current variable
    myBiomodEF <- BIOMOD_EnsembleForecasting(
      EM.output = myBiomodEM,
      projection.output = myBiomodProjTime,
      new.env = NULL,
      xy.new.env = NULL,
      selected.models = 'all',
      proj.name = NULL,
      binary.meth = NULL,
      filtered.meth = NULL,
      compress = TRUE)
  }
}

sp.names <- spData$Species
system.time(biomod_multispecies <- lapply(sp.names, MyBiomodSF))