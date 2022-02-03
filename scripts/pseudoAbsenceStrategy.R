## Pseudo-absence strategy ##
## Andre P. Silva ##
## February 2nd, 2021 ##


## Testing out different pseudo-absence strategies
myBiomodData_PArandom <- BIOMOD_FormatingData(resp.var = myResp,
                                              expl.var = env_vars,
                                              resp.xy = myRespCoord,
                                              resp.name = myRespName,
                                              PA.nb.absences = 1000,
                                              PA.strategy = c('random'),
                                              na.rm = TRUE)

myBiomodModelOut_PArandom <- BIOMOD_Modeling(myBiomodData_PArandom ,
                                             models = c('GLM'),
                                             models.options = NULL,
                                             NbRunEval = 2,
                                             DataSplit = 70,
                                             Prevalence = 0.5,
                                             VarImport = 3,
                                             models.eval.meth = c('TSS','ROC'))

myBiomodData_PAdisk <- BIOMOD_FormatingData(resp.var = myResp,
                                            expl.var = env_vars,
                                            resp.xy = myRespCoord,
                                            resp.name = myRespName,
                                            PA.nb.absences = 1000,
                                            PA.strategy = c('disk'), # collects pseudo-absences within min and max distance from occurrences
                                            PA.dist.min = 3000,
                                            PA.dist.max = 10000,
                                            na.rm = TRUE)

myBiomodModelOut_PAdisk <- BIOMOD_Modeling(myBiomodData_PAdisk,
                                           models = c('GLM'),
                                           models.options = NULL,
                                           NbRunEval = 2,
                                           DataSplit = 70,
                                           Prevalence = 0.5,
                                           VarImport = 3,
                                           models.eval.meth = c('TSS','ROC'))

get_evaluations(myBiomodModelOut_PArandom)
get_evaluations(myBiomodModelOut_PAdisk)

