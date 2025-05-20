.hiddenCall <- function(command){
  eval.parent(parse(text = command))
}

mlx.getLixoftConnectorsState <- function(quietly = TRUE) {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLixoftConnectorsState(quietly = ',quietly,')'))
  return(r)
}

mlx.initializeLixoftConnectors <- function(software = "monolix", path="", force = TRUE) {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::initializeLixoftConnectors(software = "',software ,'",
                     path = "',path,'", force=',force,')'))
  return(invisible(r))
}

mlx.setStructuralModel <- function(modelFile = NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::setStructuralModel(modelFile = "',modelFile,'")'))
}

mlx.getStructuralModel <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getStructuralModel()'))
  return(r)
}

mlx.getTests <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getTests()'))
  return(r)
}

mlx.getSAEMiterations <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSAEMiterations()'))
  return(r)
}

mlx.getPopulationParameterInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getPopulationParameterInformation()'))
  return(r)
}

mlx.getScenario <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getScenario()'))
  return(r)
}


mlx.getObservationInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getObservationInformation()'))
  return(r)
}

mlx.getLaunchedTasks <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLaunchedTasks()'))
  return(r)
}

mlx.getEstimatedLogLikelihood <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedLogLikelihood()'))
  for (k in 1:length(r)) {
    if (is.list(r[[k]]))
      r[[k]] <- unlist(r[[k]])
    if (!is.null(r[[k]]['-2LL']))
      names(r[[k]]) <- gsub("-2LL", "OFV", names(r[[k]]))
    i0 <- which(names(r[[k]])=='chosenDegree')
    if (length(i0)>0)
      r[[k]] <- r[[k]][-i0]
  }
  return(r)
}

mlx.getEstimatedPopulationParameters <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedPopulationParameters()'))
  return(r)
}

mlx.getConditionalDistributionSamplingSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getConditionalDistributionSamplingSettings()'))
  return(r)
}

mlx.getConditionalModeEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getConditionalModeEstimationSettings()'))
  return(r)
}

mlx.getContinuousObservationModel <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getContinuousObservationModel()'))
  return(r)
}

mlx.getCovariateInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getCovariateInformation()'))
  sn <- setdiff(r$name,names(r$covariate))
  if (length(sn)>0) {
    d <- mlx.getProjectSettings()$directory
    pind <- read.csv(file.path(d,"IndividualParameters/estimatedIndividualParameters.txt"))
    if (all(sn %in% names(pind)))
      r$covariate <- merge(r$covariate,pind[,c("id",sn)],by="id")
  }
  j.strat <- grep("stratification",r$type)
  if (length(j.strat) > 0) {
    strat.cov <- r$name[j.strat]
    r$covariate <- r$covariate %>% select(-strat.cov)
    r$type <- r$type[-j.strat]
    r$name <- r$name[-j.strat]
  }

  return(r)
}

mlx.getAllCovariateInformation <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getCovariateInformation()'))
  return(r)
}

mlx.getData <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getData()'))
  return(r)
}
mlx.getDemoPath <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getDemoPath()'))
  return(r)
}
mlx.getEstimatedIndividualParameters <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedIndividualParameters()'))
  return(r)
}
mlx.getEstimatedRandomEffects <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedRandomEffects()'))
  return(r)
}
mlx.getEstimatedStandardErrors <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getEstimatedStandardErrors()'))
  # for (k in 1:length(r)) {
  #   if (is.list(r[[k]]))
  #     r[[k]] <- unlist(r[[k]])
  # }
  return(r)
}
mlx.getGeneralSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getGeneralSettings()'))
  return(r)
}
mlx.getIndividualParameterModel <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getIndividualParameterModel()'))
  return(r)
}
mlx.getLogLikelihoodEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getLogLikelihoodEstimationSettings()'))
  return(r)
}
mlx.getPopulationParameterEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getPopulationParameterEstimationSettings()'))
  return(r)
}
mlx.getProjectSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getProjectSettings()'))
  return(r)
}
mlx.getSimulatedIndividualParameters <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSimulatedIndividualParameters()'))
  if (is.factor(r$rep))  r$rep <- as.numeric(as.character(r$rep))
  return(r)
}
mlx.getSimulatedRandomEffects <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getSimulatedRandomEffects()'))
  if (is.factor(r$rep))  r$rep <- as.numeric(as.character(r$rep))
  return(r)
}
mlx.getStandardErrorEstimationSettings <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getStandardErrorEstimationSettings()'))
  return(r)
}
mlx.runConditionalDistributionSampling <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::runConditionalDistributionSampling()'))
}
mlx.runConditionalModeEstimation <- function() {
  .hiddenCall(paste0('r <- lixoftConnectors::runConditionalModeEstimation()'))
}
mlx.runStandardErrorEstimation <- function(linearization=NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::runStandardErrorEstimation(linearization = ',linearization,')'))
}
mlx.runScenario <- function() {
  .hiddenCall('r <- lixoftConnectors::runScenario()')
}
mlx.setInitialEstimatesToLastEstimates <- function(fixedEffectsOnly = F) {
  .hiddenCall(paste0('r <- lixoftConnectors::setInitialEstimatesToLastEstimates(fixedEffectsOnly=fixedEffectsOnly)'))
}

mlx.setPopulationParameterInformation <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setPopulationParameterInformation(a)'))
}

mlx.loadProject <- function(projectFile=NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::loadProject(projectFile = "',projectFile,'")'))
}

mlx.setScenario <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setScenario(a)'))
}
mlx.setConditionalDistributionSamplingSettings <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setConditionalDistributionSamplingSettings(a)'))
}
mlx.setConditionalModeEstimationSettings <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setConditionalModeEstimationSettings(a)'))
}
mlx.computePredictions <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::computePredictions(a)'))
}
mlx.setCovariateModel <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setCovariateModel(a)'))
}
mlx.setIndividualParameterModel <- function(a) {
  .hiddenCall(paste0('lixoftConnectors::setIndividualParameterModel(a)'))
}

mlx.getMapping <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getMapping()'))
  return(r)
}

mlx.setMapping <- function(a) {
  .hiddenCall(paste0('lixoftConnectors::setMapping(a)'))
}

mlx.setCorrelationBlocks <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setCorrelationBlocks(a)'))
}
mlx.setData <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setData(a)'))
}
mlx.setErrorModel <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setErrorModel(a)'))
}
mlx.setGeneralSettings <- function(g) {
  .hiddenCall(paste0('r <- lixoftConnectors::setGeneralSettings(g)'))
}
mlx.setLogLikelihoodEstimationSettings <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setLogLikelihoodEstimationSettings(a)'))
}
mlx.setPopulationParameterEstimationSettings <- function(g) {
  .hiddenCall(paste0('r <- lixoftConnectors::setPopulationParameterEstimationSettings(g)'))
}
mlx.newProject <- function(data = NULL, modelFile = NULL) {
  .hiddenCall(paste0('r <- lixoftConnectors::newProject(data = data, modelFile = modelFile)'))
}

mlx.setProjectSettings <- function(directory = NULL, dataandmodelnexttoproject = NULL) {
  if (!is.null(directory)) {
    .hiddenCall(paste0('r <- lixoftConnectors::setProjectSettings(directory = directory)'))
  } else {
    if (!is.null(dataandmodelnexttoproject)) {
      .hiddenCall(paste0('r <- lixoftConnectors::setProjectSettings(dataandmodelnexttoproject = dataandmodelnexttoproject)'))
    }
  }
}
mlx.setStandardErrorEstimationSettings  <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setStandardErrorEstimationSettings (a)'))
}
mlx.setObservationDistribution  <- function(a) {
  .hiddenCall(paste0('r <- lixoftConnectors::setObservationDistribution (a)'))
}


mlx.saveProject <- function(projectFile=NULL) {
  if (is.null(projectFile)) {
    .hiddenCall(paste0('r <- lixoftConnectors::saveProject()'))
  } else {
    .hiddenCall(paste0('r <- lixoftConnectors::saveProject(projectFile = projectFile)'))
  }
}

mlx.runPopulationParameterEstimation <- function(parameters=NULL) {
  r <- NULL
  if (!is.null(parameters))
    .hiddenCall(paste0('r <- lixoftConnectors::runPopulationParameterEstimation(parameters=parameters)'))
  else
    .hiddenCall(paste0('r <- lixoftConnectors::runPopulationParameterEstimation()'))
  .hiddenCall(paste0('r0 <- lixoftConnectors::runConditionalModeEstimation()'))
  return(r)
}

mlx.runLogLikelihoodEstimation <- function(linearization = FALSE) {
  .hiddenCall(paste0('r <- lixoftConnectors::runLogLikelihoodEstimation(linearization = linearization)'))
}

mlx.getLibraryModelName <- function(library) {
  .hiddenCall(paste0('r <- lixoftConnectors::getLibraryModelName(library)'))
}

mlx.saveFormattedFile <- function(path) {
  .hiddenCall(paste0('r <- lixoftConnectors::getFormatting()'))
  .hiddenCall(paste0('r$formattedFile <- path'))
  .hiddenCall(paste0('do.call(lixoftConnectors::formatData, r)'))
}

mlx.getFormatting <- function() {
  r <- NULL
  .hiddenCall(paste0('r <- lixoftConnectors::getFormatting()'))
  return(r)
}

Rsmlx.p.weight <- function(p,pw,coef){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::p.weight(p=p,pw=pw,coef=coef)'))
  return(r)
}

print_result <- function (print, summary.file, to.cat = NULL, to.print = NULL)
{
  if (file.exists(summary.file))
    sink(summary.file, append = TRUE)
  else sink(summary.file)
  if (!is.null(to.cat))
    cat(to.cat)
  if (!is.null(to.print))
    print(to.print)
  sink()
  if (print) {
    if (!is.null(to.cat))
      cat(to.cat)
    if (!is.null(to.print))
      print(to.print)
  }
}

Rsmlx.compute.criterion <- function(criterion, method.ll, weight = NULL, pen.coef = NULL){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::compute.criterion(criterion,method.ll,weight,pen.coef)'))
  return(r)
}

Rsmlx.formatLL <- function(ll,criterion,cr,is.weight,is.prior=F){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::formatLL(ll,criterion,cr,is.weight,is.prior)'))
  return(r)
}

Rsmlx.correlationTest <- function(project = NULL, n.sample = NULL, plot = FALSE){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::correlationTest(project,n.sample,plot)'))
  return(r)
}

Rsmlx.sortCov <- function(r,cov.ini){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::sortCov(r,cov.ini)'))
  return(r)
}

Rsmlx.formatCovariateModel <- function(m,cov.ini=NULL){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::formatCovariateModel(m,cov.ini)'))
  return(r)
}

Rsmlx.formatErrorModel <- function(m){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::formatErrorModel(m)'))
  return(r)
}

Rsmlx.covariateModelSelection <- function(pen.coef = NULL, weight = 1, n.full = 10, nb.model = 1,
                                          covToTransform = NULL, covFix = NULL, direction = "both",
                                          paramToUse = "all", steps = 1000, p.max = 1, sp0 = NULL,
                                          iter = 1, correlation.model = NULL, eta = NULL){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::covariateModelSelection(pen.coef,weight,n.full,nb.model,covToTransform,covFix,direction,paramToUse,steps,p.max,sp0,iter,correlation.model,eta)'))
  return(r)
}

Rsmlx.prcheck <- function(project, f = NULL, settings = NULL, model = NULL,
                          paramToUse = NULL, parameters = NULL, level = NULL, tests = NULL,
                          nboot = NULL, method = NULL){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::prcheck(project,f,settings,model,paramToUse,parameters,level,tests,nboot,method)'))
  return(r)
}

Rsmlx.buildmlx.check <- function(project, final.project, model, paramToUse, covToTest,
                                 covToTransform, center.covariate, criterion, linearization,
                                 ll, test, direction, steps, max.iter, explor.iter, seq.cov,
                                 seq.cov.iter, seq.corr, p.max, p.min, print, nb.model, prior,
                                 weight, n.full){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::buildmlx.check(project, final.project, model, paramToUse, covToTest,
  covToTransform, center.covariate, criterion, linearization,
  ll, test, direction, steps, max.iter, explor.iter, seq.cov,
  seq.cov.iter, seq.corr, p.max, p.min, print, nb.model, prior,
  weight, n.full)'))
  return(r)
}

Rsmlx.def.variable <- function(weight = NULL, prior = NULL, criterion = NULL, fix.param0 = NULL,
                               fix.param1 = NULL){
  r <- NULL
  .hiddenCall(paste0('r <- Rsmlx:::def.variable(weight,prior,criterion,fix.param0,fix.param1)'))
  return(r)
}

Rsmlx.errorModelSelection <- function(project = NULL, pen.coef = NULL, nb.model = 1, f.min = 0.001){
  r <- NULL
  .hiddenCall(' r <- Rsmlx:::errorModelSelection(project,pen.coef,nb.model,f.min)')
  return(r)
}

Rsmlx.covariate.test <- function(cov.test, covToTest, covToTransform, paramToUse){
  r <- NULL
  .hiddenCall('r <- Rsmlx:::covariate.test(cov.test, covToTest, covToTransform, paramToUse)')
  return(r)
}

Rsmlx.correlationModelSelection <- function(e0 = NULL, pen.coef = NULL, nb.model = 1, corr0 = NULL,
                                            seqmod = TRUE, prior = NULL, cor.list = NULL, weight = NULL){
  r <- NULL
  .hiddenCall('r <- Rsmlx:::correlationModelSelection(e0,pen.coef,nb.model,corr0,seqmod,prior,cor.list,weight)')
  return(r)
}

mlx.stepAIC <- function(object, scope, scale = 0, direction = c("both", "backward",
                                                                "forward"), trace = 1, keep = NULL, steps = 1000, use.start = FALSE,
                        k = 2, weight = NULL){
  r <- NULL
  .hiddenCall('r <- Rsmlx:::mlx.stepAIC(object,scope,scale,direction,trace,keep,steps,use.start,k,weight)')
  return(r)
}
