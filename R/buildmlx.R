#' buildmlx: Automatic statistical model building
#'
#' buildmlx uses SAMBA (Stochastic Approximation for Model Building Algorithm), an iterative procedure to accelerate and optimize the process of model building by identifying at each step how best to improve some of the model components. This method allows to find the optimal statistical model which minimizes some information criterion in very few steps.
#'
#' @details
#' For covariates model building, covariate selection can be achieved by stepAIC, as the original SAMBA algorithm was implemented in Rsmlx package (Prague and Lavielle, 2020  ; Mihaljevic, 2023) and by a lasso approach enhanced by stability selection (Bodinier et al., 2023).
#'
#'
#' @param project a string: the initial Monolix project
#' @param final.project a  string: the final Monolix project (default adds "_built" to the original project)
#' @param model components of the model to optimize c("residualError", "covariate", "correlation"), (default="all")
#' @param prior list of prior probabilities for each component of the model (default=NULL)
#' @param weight list of penalty weights for each component of the model (default=NULL)
#' @param coef.w1 multiplicative weight coefficient used for the first iteration only (default=0.5)
#' @param paramToUse list of parameters possibly function of covariates (default="all")
#' @param covToTest components of the covariate model that can be modified (default="all")
#' @param covToTransform list of (continuous) covariates to be log-transformed (default="none")
#' @param center.covariate TRUE/FALSE center the covariates of the final model (default=FALSE)
#' @param criterion penalization criterion to optimize c("AIC", "BIC", "BICc", gamma) (default=BICc)
#' @param linearization TRUE/FALSE whether the computation of the likelihood is based on a linearization of the model (default=FALSE)
#' @param ll TRUE/FALSE compute the observe likelihood and the criterion to optimize at each iteration
#' @param test TRUE/FALSE perform additional statistical tests for building the model (default=FALSE)
#' @param direction for stepAIC method, method for covariate search c("full", "both", "backward", "forward"), (default="full" or "both")
#' @param steps for stepAIC method, maximum number of iteration for stepAIC (default=1000)
#' @param n.full for stepAIC method, maximum number of covariates for an exhaustive comparison of all possible covariate models (default=10)
#' @param max.iter maximum number of iterations (default=20)
#' @param explor.iter number of iterations during the exploratory phase (default=2)
#' @param fError.min minimum fraction of residual variance for combined error model (default = 1e-3)
#' @param seq.cov TRUE/FALSE whether the covariate model is built before the correlation model
#' @param seq.cov.iter number of iterations before building the correlation model (only when seq.cov=F, default=0)
#' @param seq.corr TRUE/FALSE whether the correlation model is built iteratively (default=TRUE)
#' @param p.max maximum p-value used for removing non significant relationships between covariates and individual parameters (default=0.1 for stepAIC and 1 for lasso)
#' @param p.min vector of 3 minimum p-values used for testing the components of a new model (default=c(0.075, 0.05, 0.1))
#' @param print TRUE/FALSE display the results (default=TRUE)
#' @param nb.model number of models to display at each iteration (default=1)
#' @param nfolds for lasso method, number of folds (default=10)
#' @param alpha for lasso method, the elasticnet mixing parameter, between 0 and 1. `alpha=1` is the lasso penalty, `alpha=0` is the ridge penalty.
#' @param nSS for lasso method, number of resampling iterations for stability selection.
#' @param buildMethod the method used to build the covariate model (default="lasso")
#' @param FDR_thr for lasso method, upper-bounds in FDP of calibrated stability selection (default=0.1)
#'
#' @returns a new Monolix project with a new statistical model.
#' @export
#'
#' @references Prague M, Lavielle M.  SAMBA: A novel method for fast automatic model building in nonlinear mixed-effects models.  CPT Pharmacometrics Syst Pharmacol. 2022; 11: 161-172. doi:10.1002/psp4.12742
#'
#' Bodinier B, Filippi S, Haugdahl Nøst T, Chiquet J, Chadeau-Hyam M. Automated calibration for stability selection in penalised regression and graphical models. Journal of the Royal Statistical Society Series C: Applied Statistics. 2023 ; 72: 1375–1393. doi:10.1093/jrsssc/qlad058
#'
#' Mihaljevic F (2023).  Rsmlx: R Speaks 'Monolix'. R package version2023.1.5, <https://CRAN.R-project.org/package=Rsmlx>.
#'
#' @examples
#' \dontrun{
#' project <- getMLXdir()
#'
#' res = buildmlx(project = project,
#'                buildMethod = "lasso",
#'                model='covariate',
#'                test=FALSE)
#'
#' getIndividualParameterModel()
#' }
buildmlx <- function (project = NULL,
                      final.project = NULL,
                      model = "all",
                      prior = NULL,
                      weight = NULL,
                      coef.w1 = 0.5,
                      paramToUse = "all",
                      covToTest = "all",
                      covToTransform = "none",
                      center.covariate = FALSE,
                      criterion = "BICc",
                      linearization = FALSE,
                      ll = TRUE,
                      test = FALSE,
                      direction = NULL,
                      steps = 1000,
                      n.full = 10,
                      max.iter = 20,
                      explor.iter = 2,
                      fError.min = 0.001,
                      seq.cov = FALSE,
                      seq.cov.iter = 0,
                      seq.corr = TRUE,
                      p.max = if(buildMethod=="stepAIC"){0.1}else{1},
                      p.min = c(0.075, 0.05, 0.1),
                      print = TRUE,
                      nb.model = 1,
                      nfolds=5,
                      alpha=1,
                      nSS=1000,
                      buildMethod="lasso",
                      FDR_thr=0.1)
{
  doParallel::registerDoParallel(cluster <- parallel::makeCluster(parallel::detectCores()))
  ##########################################

  ptm <- proc.time()
  initRsmlx()
  dashed.line <- "--------------------------------------------------\n"
  plain.line <- "__________________________________________________\n"
  dashed.short <- "-----------------------\n"
  plain.short <- "_______________________\n"
  op.original <- options()
  op.new <- options()
  on.exit(options(op.original))
  op.new$lixoft_notificationOptions$warnings <- 1
  options(op.new)
  RsmlxDemo1.project <- RsmlxDemo2.project <- warfarin.data <- resMonolix <- NULL
  pi <- 4 * atan(1)
  if (!is.null(project)) {
    r <- Rsmlx.prcheck(project, f = "build", paramToUse = paramToUse,
                         model = model)
    if (r$demo)
      return(r$res)
    project <- r$project
  }
  else {
    project <- mlx.getProjectSettings()$project
  }
  method.ll <- iop.ll <- pen.coef <- NULL
  r <- Rsmlx.buildmlx.check(project, final.project, model, paramToUse,
                              covToTest, covToTransform, center.covariate, criterion,
                              linearization, ll, test, direction, steps, max.iter,
                              explor.iter, seq.cov, seq.cov.iter, seq.corr, p.max,
                              p.min, print, nb.model, prior, weight, n.full)
  if (!is.null(r$change))
    return(list(change = FALSE))
  for (j in 1:length(r)) eval(parse(text = paste0(names(r)[j],
                                                  "= r[[j]]")))
  r <- Rsmlx.def.variable(weight = weight, prior = prior, criterion = criterion)
  for (j in 1:length(r)) eval(parse(text = paste0(names(r)[j],
                                                  "= r[[j]]")))
  is.weight <- weight$is.weight
  is.prior <- NULL
  final.dir <- sub(pattern = "(.*)\\..*$", replacement = "\\1",
                   final.project)
  if (dir.exists(final.dir))
    unlink(final.dir, recursive = TRUE)
  project.dir <- mlx.getProjectSettings()$directory
  if (!dir.exists(project.dir))
    dir.create(project.dir)
  buildmlx.dir <- file.path(mlx.getProjectSettings()$directory,
                            "buildmlx")
  Sys.sleep(0.1)
  if (dir.exists(buildmlx.dir))
    unlink(buildmlx.dir, recursive = TRUE)
  Sys.sleep(0.1)
  dir.create(buildmlx.dir)
  summary.file = file.path(buildmlx.dir, "summary.txt")
  Sys.sleep(0.1)
  if (!dir.exists(final.dir))
    dir.create(final.dir, recursive = TRUE)
  to.cat <- paste0("\n", dashed.line, "\nBuilding:\n")
  if (model$covariate)
    to.cat <- c(to.cat, "  -  The covariate model\n")
  if (model$correlation)
    to.cat <- c(to.cat, "  -  The correlation model\n")
  if (model$residualError)
    to.cat <- c(to.cat, "  -  The residual error model\n")
  to.cat <- c(to.cat, "\n")
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  print.line <- FALSE
  p.ini <- mlx.getPopulationParameterInformation()
  rownames(p.ini) <- p.ini$name
  ind.omega <- grep("omega_", p.ini[["name"]])
  omega <- p.ini$name[ind.omega]
  omega.ini <- p.ini[ind.omega, ]
  error.model <- mlx.getContinuousObservationModel()$errorModel
  obs.dist <- mlx.getContinuousObservationModel()$distribution
  covariate.model <- mlx.getIndividualParameterModel()$covariateModel
  cov.ini <- names(covariate.model[[1]])
  correlation.model <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id,
                              sort)
  if (length(correlation.model) == 0)
    correlation.model <- NULL
  error.model.ini <- error.model
  covariate.model.ini <- covariate.model
  correlation.model.ini <- correlation.model
  to.cat <- ("- - - Initialization - - -\n")
  if (!print.line)
    to.cat <- paste0(plain.line, to.cat)
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  if (model$covariate) {
    to.cat <- "\nCovariate model:\n"
    to.print <- Rsmlx.formatCovariateModel(covariate.model, cov.ini)
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (model$correlation) {
    to.cat <- "\nCorrelation model:\n"
    to.print <- ifelse(!is.null(correlation.model), correlation.model,
                       "NULL")
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (model$residualError) {
    to.cat <- "\nResidual error model:\n"
    to.print <- Rsmlx.formatErrorModel(error.model)
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (!mlx.getLaunchedTasks()[["populationParameterEstimation"]]) {
    to.cat <- "\nEstimation of the population parameters using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    mlx.runPopulationParameterEstimation()
  }
  if (!mlx.getLaunchedTasks()[["conditionalDistributionSampling"]]) {
    to.cat <- "Sampling of the conditional distribution using the initial model ... \n"
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    mlx.runConditionalDistributionSampling()
  }
  gi <- mlx.getSimulatedIndividualParameters()
  gi <- gi %>% filter(rep == gi$rep[nrow(gi)]) %>% select(-rep)
  lin.ll <- method.ll == "linearization"
  launched.tasks <- mlx.getLaunchedTasks()
  if (iop.ll) {
    if (!(method.ll %in% launched.tasks[["logLikelihoodEstimation"]])) {
      if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
        mlx.runConditionalModeEstimation()
      to.cat <- "Estimation of the log-likelihood of the initial model ... \n"
      print_result(print, summary.file, to.cat = to.cat,
                           to.print = NULL)
      mlx.runLogLikelihoodEstimation(linearization = lin.ll)
    }
    ll.ini <- Rsmlx.compute.criterion(criterion, method.ll, weight,
                                pen.coef)
    ll <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                   criterion, ll.ini, is.weight, is.prior)
    list.criterion <- ll.ini
    to.cat <- paste0("\nEstimated criteria (", method.ll,
                     "):\n")
    to.print <- round(ll, 2)
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  mlx.saveProject(final.project)
  if (!is.null(covToTransform)) {
    covariate <- mlx.getCovariateInformation()$covariate
    for (cov.name in covToTransform) {
      covkj <- covariate[[cov.name]]
      lckj <- paste0("logt", toupper(substr(cov.name,
                                            1, 1)), substr(cov.name, 2, nchar(cov.name)))
      if (!(lckj %in% mlx.getCovariateInformation()$name)) {
        tr.str <- paste0(lckj, " = \"log(", cov.name,
                         "/", signif(mean(covkj), digits = 2), ")\"")
        trs <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",
                      tr.str, ")")
        eval(parse(text = trs))
        foo <- colnames(weight$covariate)
        weight$covariate <- cbind(weight$covariate,
                                  weight$covariate[, cov.name])
        colnames(weight$covariate) <- c(foo, lckj)
      }
      covToTest <- c(covToTest, lckj)
    }
    covToTest <- unique(setdiff(covToTest, covToTransform))
    covToTransform <- NULL
    mlx.saveProject(final.project)
    if (!mlx.getLaunchedTasks()[["populationParameterEstimation"]]) {
      to.cat <- "\nEstimation of the population parameters using the transformed covariates ... \n"
      print_result(print, summary.file, to.cat = to.cat,
                           to.print = NULL)
      mlx.runPopulationParameterEstimation()
    }
    if (!mlx.getLaunchedTasks()[["conditionalDistributionSampling"]]) {
      to.cat <- "Sampling of the conditional distribution using the the transformed covariates ... \n"
      print_result(print, summary.file, to.cat = to.cat,
                           to.print = NULL)
      mlx.runConditionalDistributionSampling()
    }
  }
  stop.test <- FALSE
  corr.test <- FALSE
  iter <- 0
  cov.names0 <- cov.names <- NULL
  if (identical(covToTest, "all"))
    covFix = NULL
  else covFix <- setdiff(mlx.getCovariateInformation()$name,
                         covToTest)
  if (iop.ll) {
    ll.prev <- Inf
    ll.new <- ll.ini
  }
  sp0 <- NULL
  cov.test <- NULL
  e <- NULL
  while (!stop.test) {
    iter <- iter + 1
    if (iter == 1) {
      weight.covariate <- weight$covariate * coef.w1
      weight.correlation <- weight$correlation * coef.w1
    }
    else {
      weight.covariate <- weight$covariate
      weight.correlation <- weight$correlation
    }
    to.cat <- paste0(plain.line, "- - - Iteration ", iter,
                     " - - -\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    obs.dist0 <- obs.dist
    error.model0 <- error.model
    covariate.model0 <- covariate.model
    correlation.model0 <- correlation.model
    if (iop.ll)
      ll0 <- ll
    if (model$residualError) {
      res.error <- Rsmlx.errorModelSelection(pen.coef = pen.coef[-c(1,
                                                                      2)], nb.model = nb.model, f.min = fError.min)
      if (nb.model == 1)
        error.model <- res.error
      else {
        error.model <- mlx.getContinuousObservationModel()$errorModel
        for (k in (1:length(error.model))) error.model[[k]] <- as.character(res.error[[k]]$error.model[1])
      }
    }
    pmax.cov <- p.max
    if (model$covariate) {
      res.covariate <- covariateModelSelection(buildMethod = buildMethod,
                                               nfolds = nfolds, alpha = alpha, nSS = nSS,
                                               pen.coef = pen.coef[1],
                                               nb.model = nb.model, weight = weight.covariate,
                                               covFix = covFix,
                                               direction = direction, steps = steps, p.max = pmax.cov,
                                               paramToUse = paramToUse, sp0 = sp0, iter = iter,
                                               correlation.model = correlation.model, n.full = n.full,
                                               eta = e,
                                               covariate.model = covariate.model, criterion = criterion,
                                               FDR_thr = FDR_thr,print=print)
      res.covariate$res <- Rsmlx.sortCov(res.covariate$res,
                                           cov.ini)
      if (iter > explor.iter)
        sp0 <- res.covariate$sp
      covToTransform <- setdiff(covToTransform, res.covariate$tr0)
      covariate.model <- res.covariate$model
      e <- res.covariate$residuals
      if (nb.model == 1)
        cov.select <- rownames(res.covariate$res)
      else cov.select <- rownames(res.covariate$res[[1]])
      cov.names <- lapply(covariate.model[cov.select],
                          function(x) {
                            sort(names(which(x)))
                          })
      cov.names0 <- lapply(covariate.model0[cov.select],
                           function(x) {
                             sort(names(which(x)))
                           })
      cov.test <- Rsmlx.covariate.test(cov.test, covToTest,
                                         covToTransform, paramToUse)
    }
    else {
      e <- mlx.getSimulatedRandomEffects()
    }
    if (model$correlation & !corr.test) {
      if (isTRUE(all.equal(cov.names0, cov.names)))
        corr.test <- TRUE
      if (!seq.cov & iter > seq.cov.iter)
        corr.test <- TRUE
      if (corr.test & (seq.cov == TRUE | seq.cov.iter > 0)) {
        to.cat <- "Start building correlation model too ... \n"
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
      }
    }
    if (model$correlation & corr.test) {
      pen.corr <- ifelse(iter <= 1, pen.coef[1], pen.coef[1])
      res.correlation <- Rsmlx.correlationModelSelection(e0 = e,
                                                           pen.coef = pen.corr, nb.model = nb.model, corr0 = correlation.model0,
                                                           seqmod = seq.corr, weight = weight.correlation)
      if (nb.model == 1)
        correlation.model <- res.correlation
      else correlation.model <- res.correlation[[1]]$block
      if (length(correlation.model) == 0)
        correlation.model <- NULL
    }
    else {
      res.correlation <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id,
                                sort)
    }
    if (length(res.correlation) == 0)
      res.correlation <- NULL
    if (max.iter > 0 | nb.model > 1) {
      eq.cov <- isTRUE(all.equal(cov.names0, cov.names)) &
        (pmax.cov == p.max)
      eq.err <- isTRUE(all.equal(error.model0, error.model))
      eq.corr <- isTRUE(all.equal(correlation.model0,
                                  correlation.model))
      eq.dist <- isTRUE(all.equal(obs.dist0, obs.dist))
      if (!model$correlation | corr.test) {
        if (eq.cov & eq.err & eq.dist & eq.corr)
          stop.test <- TRUE
        if (model$covariate & eq.cov & eq.dist & eq.corr)
          stop.test <- TRUE
      }
      if (stop.test) {
        to.cat <- "\nNo difference between two successive iterations\n"
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
      }
      if (!stop.test | nb.model > 1) {
        if (model$covariate) {
          to.cat <- "\nCovariate model:\n"
          to.print <- res.covariate$res
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = to.print)
        }
        if (model$correlation) {
          to.cat <- "\nCorrelation model:\n"
          if (!is.null(res.correlation))
            to.print <- res.correlation
          else to.print <- "NULL"
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = to.print)
        }
        if (model$residualError) {
          to.cat <- "\nResidual error model:\n"
          to.print <- res.error
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = to.print)
        }
      }
    }
    if (!stop.test) {
      p.est <- mlx.getEstimatedPopulationParameters()
      mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = TRUE)
      p.ini <- mlx.getPopulationParameterInformation()
      rownames(p.ini) <- p.ini$name
      i.omega <- which(grepl("omega_", p.ini$name) & (!identical(p.ini$method,
                                                                 "FIXED")))
      p.ini$initialValue[i.omega] <- p.est[p.ini$name[i.omega]] *
        3
      jcor <- grep("corr_", p.ini$name)
      if (length(jcor) > 0)
        p.ini <- p.ini[-jcor, ]
      mlx.setPopulationParameterInformation(p.ini)
      if (model$residualError) {
        emodel <- error.model
        odist <- mlx.getContinuousObservationModel()$distribution
        for (k in (1:length(emodel))) {
          if (identical(emodel[[k]], "exponential")) {
            emodel[[k]] <- "constant"
            odist[[k]] <- "lognormal"
          }
          else {
            odist[[k]] <- "normal"
          }
        }
        mlx.setErrorModel(emodel)
        mlx.setObservationDistribution(odist)
      }
      if (model$covariate) {
        if (length(res.covariate$add.covariate) > 0) {
          for (k in 1:length(res.covariate$add.covariate)) eval(parse(text = res.covariate$add.covariate[[k]]))
        }
        mlx.setCovariateModel(covariate.model)
      }
      if (model$correlation & corr.test)
        mlx.setCorrelationBlocks(correlation.model)
      if (max.iter > 0) {
        if (iop.ll) {
          if (ll.new > ll.prev) {
            g <- mlx.getGeneralSettings()
            g$autochains <- FALSE
            g$nbchains <- g$nbchains + 1
            mlx.setGeneralSettings(g)
          }
          ll.prev <- ll.new
        }
        g = mlx.getConditionalDistributionSamplingSettings()
        g$nbminiterations <- max(100, g$nbminiterations)
        mlx.setConditionalDistributionSamplingSettings(g)
      }
      mlx.saveProject(final.project)
      if (dir.exists(final.dir))
        unlink(final.dir, recursive = TRUE)
      mlx.loadProject(final.project)
      if (max.iter > 0) {
        to.cat <- paste0("\nRun scenario for model ",
                         iter, " ... \nEstimation of the population parameters... \n")
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
        mlx.runPopulationParameterEstimation()
        to.cat <- "Sampling from the conditional distribution... \n"
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
        mlx.runConditionalDistributionSampling()
        gi <- mlx.getSimulatedIndividualParameters()
        gi <- gi %>% filter(rep == gi$rep[nrow(gi)]) %>%
          select(-rep)
        if (iop.ll) {
          to.cat <- "Estimation of the log-likelihood... \n"
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
          if (lin.ll & !launched.tasks[["conditionalModeEstimation"]])
            mlx.runConditionalModeEstimation()
          mlx.runLogLikelihoodEstimation(linearization = lin.ll)
          ll.new <- Rsmlx.compute.criterion(criterion, method.ll,
                                              weight, pen.coef)
          list.criterion <- c(list.criterion, ll.new)
        }
        buildmlx.project.iter <- file.path(buildmlx.dir,
                                           paste0("iteration", iter, ".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        if (iop.ll) {
          if (stop.test)
            ll <- ll0
          else {
            ll <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                                   criterion, ll.new, is.weight)
          }
          to.cat <- paste0("\nEstimated criteria (",
                           method.ll, "):\n")
          to.print <- round(ll, 2)
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = to.print)
        }
        if (iter >= max.iter) {
          stop.test <- TRUE
          to.cat <- "Maximum number of iterations reached\n"
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
        }
      }
    }
    if (max.iter == 0)
      stop.test <- TRUE
    mlx.saveProject(final.project)
  }
  change.error.model <- NULL
  if (iop.ll) {
    ll.min <- min(list.criterion)
    iter.opt <- which.min(list.criterion)
    if (iter.opt == 1)
      buildmlx.project.iter <- project
    else {
      buildmlx.project.iter <- file.path(buildmlx.dir,
                                         paste0("iteration", iter.opt - 1, ".mlxtran"))
    }
    mlx.loadProject(buildmlx.project.iter)
    mlx.saveProject(final.project)
  }
  else {
    iter.opt <- iter
  }
  if (test) {
    if (!iop.ll) {
      g <- as.list(mlx.getLaunchedTasks())$logLikelihoodEstimation
      if (!linearization & !("importanceSampling" %in%
                             g))
        mlx.runLogLikelihoodEstimation(linearization = FALSE)
      if (linearization & !("linearization" %in% g))
        mlx.runLogLikelihoodEstimation(linearization = TRUE)
      ll.min <- Rsmlx.compute.criterion(criterion, method.ll,
                                          weight, pen.coef)
    }
    to.cat <- paste0(plain.line, "- - - Further tests - - -\n")
    print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
    if (model$covariate & sum(mlx.getIndividualParameterModel()$variability$id) >
        1) {
      g0 <- mlx.getIndividualParameterModel()
      covariate <- random.effect <- p.ttest <- p.lrt <- in.model <- p.value <- NULL
      r.test <- Rsmlx.covariate.test(cov.test, covToTest, covToTransform,
                               paramToUse)
      r.test <- r.test %>% filter(!in.model) %>% select(-in.model)
      if (is.weight) {
        w.cov <- weight$covariate[cbind(r.test[["parameter"]],
                                        r.test[["covariate"]])]
        r.test <- r.test %>% mutate(p.value = Rsmlx.p.weight(p.value,
                                                               w.cov, pen.coef[1]))
      }
      r.cov0 <- res.covariate$r.cov0
      for (j in 1:nrow(r.test)) {
        pj <- r.test$parameter[j]
        if (!is.null(r.cov0[[pj]])) {
          if (r.test$covariate[j] %in% r.cov0[[pj]])
            r.test$p.value[j] <- 1
        }
      }
      i.min <- which(as.numeric(r.test$p.value) < p.min[1])
      if (length(i.min) > 0) {
        to.cat <- paste0(plain.short, "Add parameters/covariates relationships:\n")
        to.print <- r.test[i.min, ]
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = to.print)
        g <- mlx.getIndividualParameterModel()
        stop.test <- FALSE
        for (i in i.min) {
          param.i <- r.test$parameter[i]
          cov.i <- r.test$covariate[i]
          g$covariateModel[[param.i]][cov.i] <- TRUE
        }
        mlx.setIndividualParameterModel(g)
        iter <- iter + 1
        to.cat <- paste0("\nRun scenario for model ",
                         iter, " ... \nEstimation of the population parameters... \n")
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
        buildmlx.project.iter <- file.path(buildmlx.dir,
                                           paste0("iteration", iter, ".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        mlx.runPopulationParameterEstimation()
      }
      else {
        stop.test <- FALSE
      }
      if (any(mlx.getObservationInformation()$type !=
              "continuous"))
        mlx.runStandardErrorEstimation(linearization = FALSE)
      else mlx.runStandardErrorEstimation(linearization = TRUE)
      g <- mlx.getIndividualParameterModel()
      n.param <- g$name
      n.cov <- names(g$covariateModel[[1]])
      r.test <- mlx.getTests()$wald
      pv <- as.numeric(gsub("<", "", r.test$p.value))
      pv[which(is.nan(pv))] <- 0.99999
      list.ipc <- NULL
      for (np in n.param) {
        gp <- g$covariateModel[[np]]
        ngp <- names(which(gp))
        if (length(ngp) > 0) {
          for (nc in ngp) {
            g$covariateModel[[np]][nc] <- FALSE
            ipc <- grep(paste0("beta_", np, "_", nc),
                        r.test$parameter)
            pv[ipc] <- Rsmlx.p.weight(pv[ipc], weight$covariate[np,
                                                                  nc], pen.coef[1])
            if (min(pv[ipc]) < p.min[2])
              g$covariateModel[[np]][nc] <- TRUE
            else list.ipc <- c(list.ipc, ipc)
          }
        }
      }
      if (identical(g$covariateModel, g0$covariateModel))
        stop.test <- TRUE
      if (length(list.ipc) > 0) {
        to.cat <- paste0(plain.short, "Remove parameters/covariates relationships:\n")
        method <- statistics <- parameter <- NULL
        to.print <- (r.test %>% select(-c(method, statistics)) %>%
                       rename(coefficient = parameter))[list.ipc,
                       ]
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = to.print)
      }
      g1 <- g
      if (!stop.test) {
        if (length(list.ipc) > 0) {
          mlx.setIndividualParameterModel(g)
          iter <- iter + 1
          to.cat <- paste0("\nRun scenario for model ",
                           iter, " ... \nEstimation of the population parameters... \n")
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
          buildmlx.project.iter <- file.path(buildmlx.dir,
                                             paste0("iteration", iter, ".mlxtran"))
          mlx.saveProject(buildmlx.project.iter)
          mlx.runPopulationParameterEstimation()
        }
        if (lin.ll) {
          if (!launched.tasks[["conditionalModeEstimation"]])
            mlx.runConditionalModeEstimation()
        }
        else {
          to.cat <- "Sampling from the conditional distribution... \n"
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
          mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
        mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- Rsmlx.compute.criterion(criterion, method.ll,
                                            weight, pen.coef)
        ll <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                       criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (", method.ll,
                         "):\n")
        to.print <- round(ll, 2)
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = to.print)
        if (ll.new < ll.min) {
          ll.min <- ll.new
          mlx.saveProject(final.project)
        }
        else {
          mlx.loadProject(final.project)
        }
      }
      g <- as.list(mlx.getLaunchedTasks())$standardErrorEstimation
      if (!linearization & !("stochasticApproximation" %in%
                             g))
        mlx.runStandardErrorEstimation(linearization = FALSE)
      if (!("linearization" %in% g))
        mlx.runStandardErrorEstimation(linearization = TRUE)
      r.test <- mlx.getTests()$wald
      g <- mlx.getIndividualParameterModel()
      n.param <- g$name
      n.cov <- names(g$covariateModel[[1]])
      pv <- as.numeric(gsub("<", "", r.test$p.value))
      pv[which(is.nan(pv))] <- 0
      list.ipc <- NULL
      for (np in n.param) {
        gp <- g$covariateModel[[np]]
        ngp <- names(which(gp))
        if (length(ngp) > 0) {
          for (nc in ngp) {
            ipc <- grep(paste0("beta_", np, "_", nc),
                        r.test$parameter)
            pv[ipc] <- Rsmlx.p.weight(pv[ipc], weight$covariate[np,
                                                          nc], pen.coef[1])
            if (max(pv[ipc]) > p.min[2]) {
              g$covariateModel[[np]][nc] <- FALSE
              list.ipc <- c(list.ipc, ipc)
            }
          }
        }
      }
      if (length(list.ipc) > 0 & !identical(g$covariateModel,
                                            g0$covariateModel) & !identical(g$covariateModel,
                                                                            g1$covariateModel)) {
        to.cat <- paste0(plain.short, "Remove parameters/covariates relationships:\n")
        method <- statistics <- parameter <- NULL
        to.print <- (r.test %>% select(-c(method, statistics)) %>%
                       rename(coefficient = parameter))[list.ipc,
                       ]
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = to.print)
        mlx.setIndividualParameterModel(g)
        iter <- iter + 1
        to.cat <- paste0("\nRun scenario for model ",
                         iter, " ... \nEstimation of the population parameters... \n")
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
        buildmlx.project.iter <- file.path(buildmlx.dir,
                                           paste0("iteration", iter, ".mlxtran"))
        mlx.saveProject(buildmlx.project.iter)
        mlx.runPopulationParameterEstimation()
        if (lin.ll) {
          if (!launched.tasks[["conditionalModeEstimation"]])
            mlx.runConditionalModeEstimation()
        }
        else {
          to.cat <- "Sampling from the conditional distribution... \n"
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
          mlx.runConditionalDistributionSampling()
        }
        to.cat <- "Estimation of the log-likelihood... \n"
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = NULL)
        mlx.runLogLikelihoodEstimation(linearization = lin.ll)
        ll.new <- Rsmlx.compute.criterion(criterion, method.ll,
                                            weight, pen.coef)
        ll <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                               criterion, ll.new, is.weight, is.prior)
        to.cat <- paste0("\nEstimated criteria (", method.ll,
                         "):\n")
        to.print <- round(ll, 2)
        print_result(print, summary.file, to.cat = to.cat,
                             to.print = to.print)
        if (ll.new < ll.min) {
          ll.min <- ll.new
          mlx.saveProject(final.project)
        }
      }
    }
    if (model$correlation) {
      test.cor <- TRUE
      cor.block0 <- mlx.getIndividualParameterModel()$correlationBlocks$id
      while (test.cor) {
        mlx.loadProject(final.project)
        p.cortest <- NULL
        if (!mlx.getLaunchedTasks()$conditionalDistributionSampling)
          mlx.runConditionalDistributionSampling()
        r.test <- Rsmlx.correlationTest()$p.value %>% filter(!in.model) %>%
          rename(p.value = p.cortest)
        param1 <- gsub("eta_", "", r.test$randomEffect.1)
        param2 <- gsub("eta_", "", r.test$randomEffect.2)
        w.cor <- weight$correlation[cbind(param1, param2)] +
          weight$correlation[cbind(param2, param1)]
        r.test <- r.test %>% mutate(p.value = Rsmlx.p.weight(p.value,
                                                               w.cor, pen.coef[1]))
        i.min <- which(as.numeric(r.test$p.value) <
                         p.min[3])
        g <- mlx.getIndividualParameterModel()
        param.list <- unlist(g$correlationBlocks$id)
        if (length(i.min) > 0) {
          i.min <- i.min[which.min(r.test$p.value[i.min])]
          param1 <- gsub("eta_", "", r.test$randomEffect.1[i.min])
          param2 <- gsub("eta_", "", r.test$randomEffect.2[i.min])
          test.cor <- FALSE
          if (!(param1 %in% param.list) & !(param2 %in%
                                            param.list)) {
            l.block <- length(g$correlationBlocks$id) +
              1
            g$correlationBlocks$id[[l.block]] <- c(param1,
                                                   param2)
            test.cor <- TRUE
          }
          if (!(param1 %in% param.list) & (param2 %in%
                                           param.list)) {
            l.block <- grep(param2, g$correlationBlocks$id)
            g$correlationBlocks$id[[l.block]] <- c(g$correlationBlocks$id[[l.block]],
                                                   param1)
            test.cor <- TRUE
          }
          if ((param1 %in% param.list) & !(param2 %in%
                                           param.list)) {
            l.block <- grep(param1, g$correlationBlocks$id)
            g$correlationBlocks$id[[l.block]] <- c(g$correlationBlocks$id[[l.block]],
                                                   param2)
            test.cor <- TRUE
          }
          if (test.cor) {
            to.cat <- paste0(plain.short, "Add correlation:\n")
            to.print <- r.test[i.min, ]
            print_result(print, summary.file, to.cat = to.cat,
                                 to.print = to.print)
            to.print <- g$correlationBlocks$id
            print_result(print, summary.file, to.cat = NULL,
                                 to.print = to.print)
            gi <- mlx.getPopulationParameterInformation()
            gi$initialValue[which(gi$name == paste0("omega_",
                                                    param1))] <- 3 * gi$initialValue[which(gi$name ==
                                                                                             paste0("omega_", param1))]
            gi$initialValue[which(gi$name == paste0("omega_",
                                                    param2))] <- 3 * gi$initialValue[which(gi$name ==
                                                                                             paste0("omega_", param2))]
            mlx.setPopulationParameterInformation(gi)
            mlx.setIndividualParameterModel(g)
            iter <- iter + 1
            buildmlx.project.iter <- file.path(buildmlx.dir,
                                               paste0("iteration", iter, ".mlxtran"))
            mlx.saveProject(buildmlx.project.iter)
            to.cat <- paste0("\nRun scenario for model ",
                             iter, "  ... \nEstimation of the population parameters... \n")
            print_result(print, summary.file, to.cat = to.cat,
                                 to.print = NULL)
            mlx.runPopulationParameterEstimation()
            if (lin.ll) {
              if (!launched.tasks[["conditionalModeEstimation"]])
                mlx.runConditionalModeEstimation()
            }
            else {
              to.cat <- "Sampling from the conditional distribution... \n"
              print_result(print, summary.file, to.cat = to.cat,
                                   to.print = NULL)
              mlx.runConditionalDistributionSampling()
            }
            to.cat <- "Estimation of the log-likelihood... \n"
            print_result(print, summary.file, to.cat = to.cat,
                                 to.print = NULL)
            mlx.runLogLikelihoodEstimation(linearization = lin.ll)
            ll.new <- Rsmlx.compute.criterion(criterion, method.ll,
                                        weight, pen.coef)
            ll.disp <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                                        criterion, ll.new, is.weight, is.prior)
            to.cat <- paste0("\nEstimated criteria (",
                             method.ll, "):\n")
            to.print <- round(ll.disp, 2)
            print_result(print, summary.file, to.cat = to.cat,
                                 to.print = to.print)
            if (ll.new < ll.min) {
              ll.min <- ll.new
              mlx.saveProject(final.project)
            }
            else {
              test.cor <- FALSE
            }
          }
          else {
            test.cor <- FALSE
          }
        }
        else {
          test.cor <- FALSE
        }
      }
      mlx.loadProject(final.project)
      p <- mlx.getEstimatedPopulationParameters()
      if (any(mlx.getObservationInformation()$type !=
              "continuous")) {
        mlx.runStandardErrorEstimation(linearization = FALSE)
        se <- mlx.getEstimatedStandardErrors()$stochasticApproximation$se
        names(se) <- mlx.getEstimatedStandardErrors()$stochasticApproximation$parameter
      }
      else {
        mlx.runStandardErrorEstimation(linearization = TRUE)
        se <- mlx.getEstimatedStandardErrors()$linearization$se
        names(se) <- mlx.getEstimatedStandardErrors()$linearization$parameter
      }
      z <- as.numeric(p)/as.numeric(se[names(p)])
      names(z) <- names(p)
      pv <- pnorm(-abs(z)) * 2
      pv.corr <- pv[grep("corr_", names(pv))]
      if (length(which(pv.corr > p.min[2])) > 0) {
        mlx.saveProject(buildmlx.project.iter)
        ind <- mlx.getEstimatedIndividualParameters()$saem
        pv.block <- strsplit(gsub("corr_", "", names(pv.corr)),
                             "_")
        ind.mod <- mlx.getIndividualParameterModel()
        cb <- ind.mod$correlationBlocks$id
        cl <- unlist(cb)
        pv.tot <- rep(0, length(cl))
        names(pv.tot) <- cl
        for (k in 1:length(pv.corr)) pv.tot[pv.block[[k]]] <- pv.tot[pv.block[[k]]] +
          log(pv.corr[k])
        param0 <- names(which.max(pv.tot))
        cb <- lapply(cb, function(x) setdiff(x, param0))
        i.cb <- which(lapply(cb, length) == 1)
        if (length(i.cb) > 0)
          cb[[i.cb]] <- NULL
        if (!identical(cor.block0, cb)) {
          ind.mod$correlationBlocks$id <- cb
          to.cat <- paste0(plain.short, "Test correlation model:\n")
          to.print <- cb
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = to.print)
          mlx.setInitialEstimatesToLastEstimates(fixedEffectsOnly = TRUE)
          mlx.setIndividualParameterModel(ind.mod)
          iter <- iter + 1
          buildmlx.project.iter <- file.path(buildmlx.dir,
                                             paste0("iteration", iter, ".mlxtran"))
          mlx.saveProject(buildmlx.project.iter)
          to.cat <- paste0("Run scenario for model ",
                           iter, "  ... \nEstimation of the population parameters... \n")
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
          mlx.runPopulationParameterEstimation(parameters = ind)
          if (lin.ll) {
            if (!launched.tasks[["conditionalModeEstimation"]])
              mlx.runConditionalModeEstimation()
          }
          else {
            to.cat <- "Sampling from the conditional distribution... \n"
            print_result(print, summary.file, to.cat = to.cat,
                                 to.print = NULL)
            mlx.runConditionalDistributionSampling()
          }
          to.cat <- "Estimation of the log-likelihood... \n"
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = NULL)
          mlx.runLogLikelihoodEstimation(linearization = lin.ll)
          ll.new <- Rsmlx.compute.criterion(criterion, method.ll,
                                              weight, pen.coef)
          ll <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                                 criterion, ll.new, is.weight, is.prior)
          to.cat <- paste0("\nEstimated criteria (",
                           method.ll, "):\n")
          to.print <- round(ll, 2)
          print_result(print, summary.file, to.cat = to.cat,
                               to.print = to.print)
          if (ll.new < ll.min) {
            ll.min <- ll.new
            mlx.saveProject(final.project)
          }
        }
      }
    }
  }
  if (model$covariate & nb.model > 1)
    res.covariate$res <- Rsmlx.sortCov(res.covariate$res[[1]],
                                         cov.ini)
  mlx.loadProject(final.project)
  if (model$covariate)
    covariate.model.print <- Rsmlx.formatCovariateModel(mlx.getIndividualParameterModel()$covariateModel)
  if (model$correlation) {
    correlation.model.print <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id,
                                      sort)
    if (length(correlation.model.print) == 0)
      correlation.model.print <- NULL
  }
  if (model$residualError)
    error.model.print <- Rsmlx.formatErrorModel(mlx.getContinuousObservationModel()$errorModel)
  if (iop.ll) {
    ll.final <- Rsmlx.compute.criterion(criterion, method.ll,
                                          weight, pen.coef)
    ll <- Rsmlx.formatLL(mlx.getEstimatedLogLikelihood()[[method.ll]],
                           criterion, ll.final, is.weight, is.prior)
  }
  to.cat <- paste0(plain.line, "\nFinal statistical model:\n")
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  if (model$covariate) {
    to.cat <- "\nCovariate model:\n"
    to.print <- covariate.model.print
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (model$correlation) {
    to.cat <- "\nCorrelation model:\n"
    if (!is.null(correlation.model.print))
      to.print <- correlation.model.print
    else to.print <- "NULL"
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (model$residualError) {
    to.cat <- "\nResidual error model:\n"
    to.print <- error.model.print
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  if (iop.ll & max.iter > 0) {
    to.cat <- paste0("\nEstimated criteria (", method.ll,
                     "):\n")
    to.print <- round(ll, 2)
    print_result(print, summary.file, to.cat = to.cat, to.print = to.print)
  }
  test.del <- FALSE
  if (model$covariate & center.covariate) {
    foo <- lapply(res.covariate$model, function(x) {
      which(x)
    })
    cov.model <- unique(unlist(lapply(foo, function(x) {
      names(x)
    })))
    cov.type <- mlx.getCovariateInformation()$type[cov.model]
    cov.cont <- names(cov.type[cov.type == "continuous"])
    for (ck in cov.cont) {
      cck <- paste0("c", ck)
      covk <- mlx.getCovariateInformation()$covariate[[ck]]
      tr.str <- paste0(cck, " = \"", ck, "-", signif(mean(covk),
                                                     digits = 2), "\"")
      tr.str <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",
                       tr.str, ")")
      eval(parse(text = tr.str))
      g = mlx.getIndividualParameterModel()$covariateModel
      cg <- lapply(g, function(x) {
        foo <- x[ck]
        x[ck] <- x[cck]
        x[cck] = foo
        return(x)
      })
      mlx.setCovariateModel(cg)
      test.del <- TRUE
      covariate.model <- cg
    }
  }
  if (test.del & dir.exists(final.dir)) {
    unlink(final.dir, recursive = TRUE)
  }
  g = mlx.getScenario()
  g$tasks[[2]] = TRUE
  mlx.setScenario(g)
  mlx.saveProject(final.project)
  mlx.loadProject(final.project)
  error.model <- mlx.getContinuousObservationModel()$errorModel
  covariate.model <- mlx.getIndividualParameterModel()$covariateModel
  correlation.model <- lapply(mlx.getIndividualParameterModel()$correlationBlocks$id,
                              sort)
  if (length(correlation.model) == 0)
    correlation.model <- NULL
  dt <- proc.time() - ptm
  res <- list(project = final.project, niter = iter.opt, time = dt["elapsed"])
  if (model$covariate)
    res <- c(res, list(covariate.model = covariate.model))
  if (model$correlation)
    res <- c(res, list(correlation.model = correlation.model))
  if (model$residualError)
    res <- c(res, list(error.model = error.model))
  to.cat <- paste0("\ntotal time: ", round(dt["elapsed"],
                                           digits = 1), "s\n", plain.line)
  print_result(print, summary.file, to.cat = to.cat, to.print = NULL)
  res$change <- !(identical(error.model, error.model.ini) &
                    identical(covariate.model, covariate.model.ini) & identical(correlation.model,
                                                                                correlation.model.ini))
  res$change.error.model <- change.error.model
  res$weight <- weight
  res$covToTest <- covToTest
  options(op.original)
  file.copy(from = summary.file, to = file.path(final.dir,
                                                "summary_buildmlx.txt"))
  return(list(res=res,time=round(dt["elapsed"], digits = 1),iter=iter))
}



# covariate Model selection -----------------------------------------------


covariateModelSelection <- function(buildMethod,
                                    nfolds = 5,
                                    alpha = 1,
                                    nSS=1000,
                                    covFix = NULL,
                                    pen.coef=NULL,
                                    weight=1,
                                    n.full = 10,
                                    nb.model = 1,
                                    direction="both",
                                    paramToUse="all",
                                    eta=NULL,
                                    p.max=1,
                                    steps=1000,
                                    sp0=NULL,
                                    iter=1,
                                    correlation.model=NULL,
                                    covariate.model=NULL,
                                    criterion = "BIC",
                                    FDR_thr=0.10,
                                    print=TRUE){
  if(buildMethod %in% c("stepAIC")){
    Rsmlx.covariateModelSelection(pen.coef = pen.coef,
                                    nb.model = nb.model, weight = weight,
                                    covToTransform = NULL, covFix = covFix,
                                    direction = direction, steps = steps, p.max = p.max,
                                    paramToUse = paramToUse, sp0 = sp0, iter = iter,
                                    correlation.model = correlation.model, n.full = n.full,
                                    eta = eta)
  }else if(buildMethod=="lasso"){
    covariateModelSelection.lasso(nfolds,alpha,covFix,pen.coef,weight,paramToUse,eta,p.max,sp0,nSS,covariate.model,criterion,iter,FDR_thr,print=print)
  }
}



# lasso -------------------------------------------------------------------


covariateModelSelection.lasso <- function(nfolds = 5,
                                          alpha = 1,
                                          covFix = NULL,
                                          pen.coef=NULL,
                                          weight=1,
                                          paramToUse="all",
                                          eta=NULL,
                                          p.max=1,
                                          sp0=NULL,
                                          nSS=1000,
                                          covariate.model=NULL,
                                          criterion="BIC",
                                          iter=1,
                                          FDR_thr=0.10,
                                          print=TRUE){
  # Simulate Individual Parameters and setup parameters
  sp.df <- mlx.getSimulatedIndividualParameters()
  if (is.null(sp.df$rep))
    sp.df$rep <- 1
  if (!is.null(sp0)) {
    nrep0 <- max(sp0$rep)
    sp.df$rep <- sp.df$rep + nrep0
    dn <- setdiff(names(sp.df), names(sp0))
    sp0[dn] <- sp.df[dn]
    sp.df <- rbind(sp0, sp.df)
  }
  nrep <- max(sp.df$rep)
  ind.dist <- mlx.getIndividualParameterModel()$distribution
  param.names <- names(ind.dist)
  n.param <- length(param.names)
  cov.info <- mlx.getCovariateInformation()
  cov.names <- cov.info$name
  cov.types <- cov.info$type
  tcov.names <- NULL
  covariates <- cov.info$covariate
  cov.cat <- cov.names[cov.types == "categorical"]
  covariates[cov.cat] <- lapply(covariates[cov.cat], as.factor)
  indvar <- mlx.getIndividualParameterModel()$variability$id
  indvar[setdiff(param.names, paramToUse)] <- FALSE

  cov.model <- mlx.getIndividualParameterModel()$covariateModel
  r <- res <- r.cov0 <- list()
  eps <- 1e-15

  # Y transformation
  Y=sp.df[,c("rep","id",names(indvar)[which(indvar==TRUE)])]

  # Update eta cov0 and pweight
  eta.list =list()
  cov0.list = list()
  pweight.list = list()
  for(j in 1:n.param){
    dj <- ind.dist[j]
    nj <- names(dj)
    if (indvar[j]) {
      if (tolower(dj) == "lognormal") {
        Y[,nj] <- log(Y[,nj] + eps)
        colnames(Y)[which(colnames(Y)==nj)] <- paste0("log.", nj)
      }else if (tolower(dj) == "logitnormal") {
        Y[,nj] <- log((Y[,nj] + eps)/(1 - Y[,nj] + eps))
        colnames(Y)[which(colnames(Y)==nj)] <- paste0("logit.", nj)
      }else if (tolower(dj) == "probitnormal") {
        Y[,nj] <- qnorm(Y[,nj])
        colnames(Y)[which(colnames(Y)==nj)]  <- paste0("probit.", nj)
      }
      if (!is.null(eta)) {
        eta.list <- append(eta.list,list(eta[paste0("eta_", nj)]))
      }else {eta.list <- append(eta.list,list(NULL))}
      names(eta.list)[length(eta.list)] <- nj
      if (length(covFix) > 0){
        cmj <- cov.model[[nj]][covFix]
        cov0.list <- append(cov0.list,list(names(which(!cmj))))
      }else {
        cov0.list <-  append(cov0.list,list(NULL))
      }
      names(cov0.list)[length(cov0.list)] <- nj
      pweight.list <- append(pweight.list,list(weight[nj, ]))
      names(pweight.list)[length(pweight.list)] <- nj
    }
  }
  cov0.list <- updateCov0(Y, eta.list, covariates, p.max, covFix,
                          pen.coef, pweight.list, cov0.list)

  if(length(param.names[which(indvar)])>1){
    Sigma=diag(mlx.getEstimatedPopulationParameters()[paste0("omega_",param.names[which(indvar)])]**2)
  }else{
    Sigma = matrix(mlx.getEstimatedPopulationParameters()[paste0("omega_",param.names[which(indvar)])]**2,ncol=1,nrow=1)
  }
  colnames(Sigma) <- param.names[which(indvar)]
  rownames(Sigma) <- param.names[which(indvar)]

  N = length(unique(Y$id))

  Y.mat = sapply(Y[,-c(1,2),drop=FALSE],function(x){rowMeans(matrix(x,nrow=N))}) #1 : rep 2 : id
  X.mat  = covariates[,setdiff(colnames(covariates),"id")]

  p = NULL
  r.var = foreach(p = names(indvar)[which(indvar)],.export = c("applyMethodLasso","modelFromSelection"),.packages = c("ggplot2","gghighlight")) %dopar% {
    aux=applyMethodLasso(Y.mat[,stringr::str_detect(colnames(Y.mat),p),drop=FALSE],
                         X.mat,
                         Sigma[p,p],
                         cov0=cov0.list[[p]],
                         p.name=p,
                         nfolds=nfolds,
                         alpha=alpha,
                         nSS=nSS,
                         criterion=criterion,
                         covariate.model = covariate.model[[p]],
                         n_cores = max(floor(parallel::detectCores()/length(names(indvar)[which(indvar)])),1),
                         iter=iter,
                         FDR_thr=FDR_thr)

    return(aux)
  }

  to.cat = unlist(sapply(r.var,FUN=function(r){r$to.cat}))
  if(print){cat(to.cat)}

  r.var <- lapply(r.var,FUN=function(r){r[-which(names(r)=="to.cat")]})

  r <- res <- r.cov0 <- list()
  for(j in 1:n.param){
    nj=param.names[j]
    if(indvar[j]){
      ind = which(unlist(lapply(r.var,FUN=function(x){x$p.name}))==nj)
      r[[j]] <- r.var[[ind]]

      res[[j]] <- r[[j]]$res
      names(res)[j] <- param.names[j]
      r.cov0 <- append(r.cov0,list(r[[j]]$cov0))
      names(r.cov0)[length(r.cov0)] <- param.names[j]

    }else{
      r[[j]] <- list(model="fixed")
      res[[j]] <- "none"
      names(res)[j] <- param.names[j]
      r[[j]]$p.name <- nj
    }
  }

  e <- as.data.frame(lapply(r[indvar], function(x) {
    x$model$residuals
  }))
  e.names <- unlist(lapply(r[indvar], function(x) {
    x$p.name
  }))
  names(e) <- paste0("eta_", e.names)
  if (!is.null(sp.df["id"]))
    e <- cbind(sp.df["id"], e)
  if (!is.null(sp.df["rep"]))
    e <- cbind(sp.df["rep"], e)

  covariate.model <- mlx.getIndividualParameterModel()$covariateModel
  covariate <- mlx.getCovariateInformation()$covariate
  js <- 0
  trs <- list()
  tr0 <- NULL
  for (k in (1:n.param)) {
    if (!identical(res[[k]], "none")) {
      covariate.model[[k]][1:length(covariate.model[[k]])] <- FALSE
      if (indvar[k]) {
        ck <- attr(r[[k]]$model$terms, "term.labels")
        if (length(ck) > 0) {
          for (j in (1:length(ck))) {
            ckj <- ck[j]
            if (identical(substr(ckj, 1, 4), "log.")) {
              js <- js + 1
              ckj.name <- sub("log.", "", ckj)
              covkj <- covariate[[ckj.name]]
              lckj <- paste0("l", ckj.name)
              tr.str <- paste0(lckj, " = \"log(", ckj.name,
                               "/", signif(mean(covkj), digits = 2),
                               ")\"")
              trs[[js]] <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",
                                  tr.str, ")")
              tr0 <- unique(c(tr0, ckj.name))
              covariate.model[[k]][lckj] <- TRUE
            }
            else {
              covariate.model[[k]][ckj] <- TRUE
            }
          }
        }
      }
    }
  }
  res <- Rsmlx.formatCovariateModel(res)

  return(list(model = covariate.model, residuals = e, res = res,
              add.covariate = trs, sp = sp.df, tr0 = tr0, r.cov0 = r.cov0))
}

applyMethodLasso <- function(Y,X,omega,cov0,
                             nfolds=5,
                             alpha=1,
                             nSS=1000,
                             criterion="BIC",
                             covariate.model=NULL,
                             p.name=NULL,
                             n_cores = 1,
                             iter=1,
                             FDR_thr=0.10){
  to.cat = c()
  lambda=NULL

  if(criterion %in% c("BIC","BICc")){
    critFUN <- BIC
  }else if(criterion=="AIC"){
    critFUN <- AIC
  }else{
    critFUN = function(mod){criterion*length(coef(mod))}
  }

  cov.names = colnames(X)
  tparam.names = colnames(Y)

  if(!is.matrix(Y)){
    Yaux <- as.matrix(Y)
  }else{Yaux=Y}
  if(!is.matrix(X)){
    Xaux <- as.matrix(X)
  }else{Xaux=X}

  Xsc <- scale(apply(Xaux,2,FUN=as.numeric))

  if(!is.null(omega)){ rootInvOmega = 1/((omega)**(1/2)) }else{ rootInvOmega = 1 }
  Ywh <- Yaux %*% rootInvOmega
  Xwh <- kronecker(t(rootInvOmega),Xsc)
  colnames(Xwh) <- cov.names

  if(!is.null(covariate.model)){
    if(any(!(cov.names %in% names(covariate.model)))){
      savedSelection <- setNames(rep(0,length(cov.names)),cov.names)
      savedSelection[names(covariate.model)] <- as.numeric(covariate.model)
    }else{
      savedSelection = setNames(as.numeric(covariate.model),names(covariate.model))
    }
    prevSelection = covariate.model

    if(all(!prevSelection)){
      oldCriterion =critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(prevSelection)[which(prevSelection)]]
      oldCriterion = critFUN(lm(Ywh ~ Xkeep))
    }
    to.cat <- c(to.cat,paste0("\n Lasso selection, calibrated using sharp method, improving the ",criterion," criterion for ",p.name," :\n "))
    to.cat <- c(to.cat,paste0("       -> Old Criterion : ",round(oldCriterion,digits=2)),"\n")
  }

  if(is.null(cov0)){
    exclude = NULL
  }else{
    exclude = which(cov.names %in% cov0)
  }

  if(!is.null(exclude) && ncol(Xwh)-length(exclude)==0){
    selection = rep(0,ncol(Xwh))

    to.cat.here = ""
  }else if(!is.null(exclude) && ncol(Xwh)-length(exclude)==1){
    selection = rep(0,ncol(Xwh))
    selection[-exclude] <- 1
    if(all(!as.logical(selection))){
      newcriterion = critFUN(lm(Ywh ~ NULL))
    }else{
      Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
      newcriterion = critFUN(lm(Ywh~Xkeep))
    }

    to.cat.here = ""
  }else{
    VariableSelection.outputs = sharp::VariableSelection(Xwh,Ywh,exclude=exclude,nfolds=nfolds,pi_list=seq(0.50,0.99,0.01),alpha=alpha,K=nSS,n_cores=n_cores,FDP_thr = FDR_thr)

    pi_list = VariableSelection.outputs$params$pi_list
    lambda_list = VariableSelection.outputs$Lambda


    Score = VariableSelection.outputs$S_2d
    argmax_id = which(!is.na(Score),arr.ind=TRUE)
    if(nrow(argmax_id)!=0){
      resSharp = lapply(split(argmax_id,1:nrow(argmax_id)),FUN=function(arg_id){
        selection = sharp::SelectedVariables(VariableSelection.outputs,argmax_id = arg_id)
        if(all(!as.logical(selection))){
          newcriterion = critFUN(lm(Ywh ~ NULL))
        }else{
          Xkeep = Xwh[,names(selection)[which(as.logical(selection))]]
          newcriterion = critFUN(lm(Ywh~Xkeep))
        }

        if(newcriterion==-Inf){
          newcriterion = oldCriterion + 1
        }

        return(list(selection=selection,criterion=newcriterion))
      })

      indMax = which.min(sapply(resSharp,FUN=function(r){r$criterion}))
      selection = resSharp[[indMax]]$selection
      newcriterion = resSharp[[indMax]]$criterion

      df = data.frame()
      for(i in 1:ncol(Score)){
        df <- rbind(df,data.frame(lambda = signif(lambda_list,digits=2),pi = pi_list[i],Score=Score[,i]))
      }

      to.cat.here =  paste0("\n              > parameter values : ",
                            paste0(c("lambda","thresholds"),"=",c(signif(lambda_list[argmax_id[indMax,1]],3),signif(pi_list[argmax_id[indMax,2]],2)),collapse=", "))
    }else{
      selection = savedSelection
      newcriterion = oldCriterion
    }
  }

  if(newcriterion >= oldCriterion){
    to.cat <- c(to.cat,paste0("        No model improving the criterion as been find, the previous covariate model is kept."))
    selection = savedSelection
  }else{
    to.cat <- c(to.cat,paste0("        -> New Criterion : ",round(newcriterion,digits=2)))
    to.cat <- c(to.cat,to.cat.here)
  }

  model.list = modelFromSelection(Y,X,selection)


  to.cat <- c(to.cat,"\n")
  return(list(model=model.list,res=selection,cov0=cov0,p.name=p.name,to.cat = to.cat))
}


# stepAIC -----------------------------------------------------------------

covariateModelSelection.reg <-
  function (pen.coef = NULL, weight = 1, n.full = 10, nb.model = 1,
            covFix = NULL, direction = "both",
            paramToUse = "all", steps = 1000, p.max = 1, sp0 = NULL,
            iter = 1, correlation.model = NULL, eta = NULL)
  {
    sp.df <- mlx.getSimulatedIndividualParameters()
    if (is.null(sp.df$rep))
      sp.df$rep <- 1
    if (!is.null(sp0)) {
      nrep0 <- max(sp0$rep)
      sp.df$rep <- sp.df$rep + nrep0
      dn <- setdiff(names(sp.df), names(sp0))
      sp0[dn] <- sp.df[dn]
      sp.df <- rbind(sp0, sp.df)
    }
    nrep <- max(sp.df$rep)
    ind.dist <- mlx.getIndividualParameterModel()$distribution
    param.names <- names(ind.dist)
    n.param <- length(param.names)
    cov.info <- mlx.getCovariateInformation()
    cov.names <- cov.info$name
    cov.types <- cov.info$type
    tcov.names <- NULL
    covariates <- cov.info$covariate
    cov.cat <- cov.names[cov.types == "categorical"]
    covariates[cov.cat] <- lapply(covariates[cov.cat], as.factor)
    indvar <- mlx.getIndividualParameterModel()$variability$id
    indvar[setdiff(param.names, paramToUse)] <- FALSE
    cov.model <- mlx.getIndividualParameterModel()$covariateModel
    r <- res <- r.cov0 <- list()
    eps <- 1e-15
    for (j in (1:n.param)) {
      dj <- ind.dist[j]
      nj <- names(dj)
      if (indvar[j]) {
        yj <- sp.df[nj]
        if (tolower(dj) == "lognormal") {
          yj <- log(yj + eps)
          names(yj) <- paste0("log.", nj)
        }
        else if (tolower(dj) == "logitnormal") {
          yj <- log((yj + eps)/(1 - yj + eps))
          names(yj) <- paste0("logit.", nj)
        }
        else if (tolower(dj) == "probitnormal") {
          yj <- qnorm(yj)
          names(yj) <- paste0("probit.", nj)
        }
        if (!is.null(eta))
          etaj <- eta[paste0("eta_", nj)]
        else etaj <- NULL
        if (length(covFix) > 0) {
          cmj <- cov.model[[nj]][covFix]
          cov0 <- names(which(!cmj))
          cov1 <- names(which(cmj))
        }
        else {
          cov0 <- cov1 <- NULL
        }
        pwj <- weight[nj, ]
        r[[j]] <- lm.all(nj, yj, etaj, covariates, tcov.names,
                         pen.coef = pen.coef, nb.model = nb.model, pw = pwj,
                         n.full = n.full, direction = direction, steps = steps,
                         p.max = p.max, cov0 = cov0, cov1 = cov1, iter = iter)
        res[[j]] <- r[[j]]$res
        r.cov0 <- append(r.cov0, list(r[[j]]$cov0))
        names(r.cov0)[length(r.cov0)] <- param.names[j]
        names(res[[j]]) <- gsub("log[.]", "l", names(res[[j]]))
      }
      else {
        r[[j]] <- list(model = "fixed")
        res[[j]] <- "none"
      }
      r[[j]]$p.name <- nj
    }
    names(res) <- param.names
    names(r.cov0) <- names(indvar)[which(indvar)]
    e <- as.data.frame(lapply(r[indvar], function(x) {
      x$model$residuals
    }))
    e.names <- unlist(lapply(r[indvar], function(x) {
      x$p.name
    }))
    names(e) <- paste0("eta_", e.names)
    if (!is.null(sp.df["id"]))
      e <- cbind(sp.df["id"], e)
    if (!is.null(sp.df["rep"]))
      e <- cbind(sp.df["rep"], e)
    if (length(correlation.model) > 0) {
      for (ic in (1:length(correlation.model))) {
        pic <- correlation.model[[ic]]
        if (all(pic %in% e.names)) {
          ipic <- match(pic, gsub("eta_", "", names(e)))
          epic <- e[, ipic]
          gic <- solve(cov(epic))
          jic <- match(pic, names(ind.dist))
          for (itk in 1:2) {
            jk <- 0
            for (j in jic) {
              jk <- jk + 1
              dj <- ind.dist[j]
              nj <- names(dj)
              yj <- sp.df[nj]
              pwj <- weight[nj, ]
              if (tolower(dj) == "lognormal") {
                yj <- log(yj)
                names(yj) <- paste0("log.", nj)
              }
              else if (tolower(dj) == "logitnormal") {
                yj <- log(yj/(1 - yj))
                names(yj) <- paste0("logit.", nj)
              }
              else if (tolower(dj) == "probitnormal") {
                yj[[1]] <- qnorm(yj[[1]])
                names(yj) <- paste0("probit.", nj)
              }
              if (!is.null(eta))
                etaj <- eta[paste0("eta_", nj)]
              else etaj <- NULL
              if (length(covFix) > 0) {
                cmj <- cov.model[[nj]][covFix]
                cov0 <- names(which(!cmj))
                cov1 <- names(which(cmj))
              }
              else {
                cov0 <- cov1 <- NULL
              }
              ejc <- as.matrix(epic[, -jk]) %*% matrix(gic[jk,
                                                           -jk], ncol = 1)/gic[jk, jk]
              yjc <- yj + ejc
              r[[j]] <- lm.all(nj, yjc, etaj, covariates,
                               tcov.names, pen.coef = pen.coef, nb.model = nb.model,
                               pw = pwj, n.full = n.full, direction = direction,
                               steps = steps, p.max = p.max, cov0 = cov0,
                               cov1 = cov1, iter = iter)
              res[[j]] <- r[[j]]$res
              r.cov0[[j]] <- r[[j]]$cov0
              r[[j]]$p.name <- nj
              e[paste0("eta_", nj)] <- epic[, jk] <- r[[j]]$model$residuals -
                ejc
            }
          }
        }
      }
    }
    covariate.model <- mlx.getIndividualParameterModel()$covariateModel
    covariate <- mlx.getCovariateInformation()$covariate
    js <- 0
    trs <- list()
    tr0 <- NULL
    for (k in (1:n.param)) {
      if (!identical(res[[k]], "none")) {
        covariate.model[[k]][1:length(covariate.model[[k]])] <- FALSE
        if (indvar[k]) {
          ck <- attr(r[[k]]$model$terms, "term.labels")
          if (length(ck) > 0) {
            for (j in (1:length(ck))) {
              ckj <- ck[j]
              if (identical(substr(ckj, 1, 4), "log.")) {
                js <- js + 1
                ckj.name <- sub("log.", "", ckj)
                covkj <- covariate[[ckj.name]]
                lckj <- paste0("l", ckj.name)
                tr.str <- paste0(lckj, " = \"log(", ckj.name,
                                 "/", signif(mean(covkj), digits = 2),
                                 ")\"")
                trs[[js]] <- paste0("lixoftConnectors::addContinuousTransformedCovariate(",
                                    tr.str, ")")
                tr0 <- unique(c(tr0, ckj.name))
                covariate.model[[k]][lckj] <- TRUE
              }
              else {
                covariate.model[[k]][ckj] <- TRUE
              }
            }
          }
        }
      }
    }
    res <- Rsmlx.formatCovariateModel(res)
    return(list(model = covariate.model, residuals = e, res = res,
                add.covariate = trs, sp = sp.df, tr0 = tr0, r.cov0 = r.cov0))
  }


lm.all <-
  function (ny, y, eta, x, tr.names = NULL, pen.coef = NULL, nb.model = NULL,
            pw = NULL, n.full = 10, direction = "both", steps = 1000,
            p.max = 1, cov0 = NULL, cov1 = NULL, iter = 1)
  {
    N <- length(unique(x$id))
    nrep <- nrow(y)/N
    if (p.max <= 1) {
      nx <- setdiff(names(x), c("id", "rep"))
      yc <- rowMeans(matrix(y[[1]], nrow = N))
      xc <- x
      lm0 <- lm(yc ~ 1)
      nxc <- setdiff(nx, cov0)
      pjc <- NULL
      for (nc in nxc) {
        lmc <- lm(yc ~ xc[[nc]])
        pc <- signif(anova(lm0, lmc)$`Pr(>F)`[2], 4)
        pjc <- c(pjc, pc)
      }
      pjc <- Rsmlx.p.weight(pjc, pw[nxc], pen.coef)
      names(pjc) <- nxc
      if (!is.null(eta)) {
        etac <- rowMeans(matrix(eta[[1]], nrow = N))
        lm0 <- lm(etac ~ 1)
        pjec <- NULL
        for (nc in nxc) {
          lmc <- lm(etac ~ xc[[nc]])
          pc <- signif(anova(lm0, lmc)$`Pr(>F)`[2], 4)
          pjec <- c(pjec, pc)
        }
        pjec <- Rsmlx.p.weight(pjec, pw[nxc], pen.coef)
        names(pjec) <- nxc
        pjc <- pmin(pjc, pjec, na.rm = TRUE)
      }
      pjc[names(which(mlx.getIndividualParameterModel()$covariateModel[[ny]]))] <- 0
      list.c <- which(pjc > p.max)
      cov0 <- c(cov0, nxc[list.c])
      direction <- ifelse(length(setdiff(nx, cov0)) <= n.full,
                          "full", direction)
    }
    else list.c <- NULL
    x$id <- x$rep <- NULL
    nx <- ncol(x)
    l <- x
    s <- rep("0:1", nx)
    names(s) <- names(x)
    j.num <- which(!sapply(x, is.factor))
    j.num <- NULL
    if (length(j.num) > 0) {
      if (length(j.num) == 1)
        j0 <- which(min(x[, j.num]) > 0)
      else j0 <- which(sapply(x[, j.num], min) > 0)
      j0.num <- j.num[j0]
      s[j0.num] <- "0:2"
      l[, j0.num] <- log(x[, j0.num])
      names(l)[j0.num] <- paste0("log.", names(x)[j0.num])
      l[, j.num] <- scale(l[j.num], scale = FALSE)
    }
    else j0.num <- vector(length = 0)
    if (direction != "full") {
      l_data = data.frame(y = y[[1]])
      j0.num = j0.num[!names(j0.num) %in% tr.names]
      l_data = cbind.data.frame(l_data, x, l[, j0.num, drop = FALSE])
      llk = tryCatch({
        f.sature <- "y ~ ."
        if (length(cov0) > 0)
          f.sature <- paste0(f.sature, "-", paste(cov0,
                                                  collapse = "-"))
        f.sature <- as.formula(f.sature)
        model.sature = lm(f.sature, l_data)
        f.cst <- "y ~ 1"
        if (length(cov1) > 0)
          f.cst <- paste0(f.cst, "+", paste(cov1, collapse = "+"))
        f.cst <- as.formula(f.cst)
        model.cst = lm(f.cst, l_data)
        if (direction == "backward") {
          lm.sat = mlx.stepAIC(model.sature, direction = "backward",
                                       trace = FALSE, k = nrep * pen.coef, scope = list(upper = model.sature,
                                                                                        lower = model.cst), steps = steps, weight = pw)
        }
        else {
          c.cur <- names(which(mlx.getIndividualParameterModel()$covariateModel[[ny]]))
          f.cur <- "y ~ 1"
          if (length(c.cur) > 0)
            f.cur <- paste0(f.cur, "+", paste(c.cur, collapse = "+"))
          f.cur <- as.formula(f.cur)
          model.cur = lm(f.cur, l_data)
          lm.cst = MASS::stepAIC(model.cur, direction = direction,
                                 trace = FALSE, k = nrep * pen.coef, scope = list(upper = model.sature,
                                                                                  lower = model.cst), steps = steps)
        }
      }, error = function(e) {
        print("Error in stepAIC")
        return(-Inf)
      })
      Gnames = names(x)
      G <- data.frame(matrix(0, ncol = length(Gnames), nrow = 1))
      colnames(G) <- Gnames
      usedcovariates = names(llk$model)[-1]
      G[1, usedcovariates] <- 1
      ll = logLik(llk)/nrep
      df = length(coef(llk)) - 1
      criterion = -2 * ll + pen.coef * sum(pw[usedcovariates])
      res <- data.frame(ll = round(ll, digits = 3), df = df,
                        criterion = round(criterion, digits = 3))
      j0.num <- j0.num[!(names(j0.num) %in% tr.names)]
      if (length(j0.num) > 0) {
        G[names(l)[j0.num]] <- 0
        i2 <- (G[names(x)[j0.num]] == 2)
        G[names(l)[j0.num]][i2] <- 1
        G[names(x)[j0.num]][i2] <- 0
      }
      res <- cbind(G == 1, res)
      if (nb.model == 1)
        res[, c("ll", "df", "criterion")] <- NULL
      else res[2:nb.model, ] <- res[1, ]
      row.names(res) <- 1:nrow(res)
      return(list(model = llk, res = res, cov0 = cov0))
    }
    if (length(tr.names) > 0) {
      j.newc <- which((names(x) %in% tr.names))
      if (length(j.newc) > 0)
        s[j.newc] <- "0:1"
    }
    if (!is.null(cov0))
      s[cov0] <- "0"
    s <- paste(s, collapse = ",")
    s <- paste0("G <- expand.grid(", s, ")")
    eval(parse(text = s))
    if (length(tr.names) > 0) {
      jc <- which(!(names(x) %in% tr.names))
      c.names <- names(x)[jc]
      for (k in 1:length(j.newc)) {
        jk <- which(paste0("l", c.names) == tr.names[k])
        if (length(jk) > 0) {
          j <- which(names(x) == c.names[jk])
          tj <- which(names(x) == tr.names[k])
          i0 <- which(G[, j] > 0 & G[, tj] > 0)
          G <- G[-i0, ]
        }
      }
    }
    names(G) <- names(x)
    if (length(cov0) > 0) {
      i0 <- which(rowSums(G[cov0]) == 0)
      G <- G[i0, ]
    }
    if (length(cov1) > 0) {
      i1 <- which(rowSums(G[cov1] == 1) == length(cov1))
      G <- G[i1, ]
    }
    ng <- nrow(G)
    d <- ncol(G)
    ll <- df <- bic <- bic.cor <- NULL
    corb <- log(iter^2/(iter^2 + 3))
    iop.mean <- FALSE
    if (iop.mean) {
      yg <- colMeans(matrix(y[[1]], nrow = nrep))
      xg <- x
      if (nrow(l) > 0)
        lg <- l[seq(1, nrow(l), by = nrep), ]
      nrepg <- 1
    }
    else {
      yg <- y[[1]]
      xg <- x[rep(1:N, nrep), ]
      lg <- l
      nrepg <- nrep
    }
    for (k in 1:ng) {
      xk <- data.frame(y = yg)
      Gk <- G[k, , drop = FALSE]
      pwk <- pw[names(Gk)]
      j1 <- which(Gk == 1)
      if (length(j1) > 0)
        xk[names(x)[j1]] <- xg[j1]
      j2 <- which(Gk == 2)
      if (length(j2) > 0)
        xk[names(l)[j2]] <- lg[j2]
      llk = tryCatch({
        lmk <- lm(y ~ ., data = xk)
        logLik(lmk)[1]/nrepg
      }, error = function(e) {
        return(-Inf)
      })
      dfk <- sum(Gk > 0)
      bick <- -2 * llk + pen.coef * sum(Gk * pwk)
      ll <- c(ll, llk)
      df <- c(df, dfk)
      bic <- c(bic, bick)
      bic.cor <- c(bic.cor, bick)
    }
    bic <- round(bic.cor, digits = 3)
    i0 <- rep(1, ng)
    mG <- ncol(G)
    for (k in seq_len(ng - 1)) {
      if (i0[k] == 1) {
        ik <- which(bic[(k):ng] == bic[k]) + k - 1
        sk <- .rowSums(G[ik, ] == 2, n = length(ik), m = mG)
        ik0 <- ik[which(sk == 0)]
        if (length(ik0) == 0)
          ik0 <- ik[order(sk)[1]]
        i0[ik] <- 0
        i0[ik0] <- 1
        i0[k] <- 1
      }
    }
    res <- data.frame(ll = round(ll, digits = 3), df = df, criterion = bic)
    res <- res[i0 == 1, ]
    G <- G[i0 == 1, , drop = FALSE]
    bic <- bic[i0 == 1]
    eval(parse(text = paste0(names(y), " <- y[[1]]")))
    obic <- order(bic)
    k.min <- obic[1]
    Gkmin <- G[k.min, ]
    j1 <- which(Gkmin == 1)
    j2 <- which(Gkmin == 2)
    if (length(j1) > 0) {
      for (k in (1:length(j1))) eval(parse(text = paste0(names(x)[j1[k]],
                                                         " <- rep(x[[j1[k]]], nrep)")))
    }
    if (length(j2) > 0) {
      for (k in (1:length(j2))) eval(parse(text = paste0(names(l)[j2[k]],
                                                         " <- l[[j2[k]]]")))
    }
    list.x <- c("1", names(x)[j1], names(l)[j2])
    form1 <- paste0(names(y), "~", paste(list.x, collapse = "+"))
    eval(parse(text = paste0("lm.min <- lm(", form1, ")")))
    lm.min$covsel = Gkmin
    nb.model0 <- min(nb.model, length(bic))
    res <- res[obic[1:nb.model0], ]
    G <- G[obic[1:nb.model], , drop = FALSE]
    j0.num <- j0.num[!(names(j0.num) %in% tr.names)]
    if (length(j0.num) > 0) {
      G[names(l)[j0.num]] <- 0
      i2 <- (G[names(x)[j0.num]] == 2)
      G[names(l)[j0.num]][i2] <- 1
      G[names(x)[j0.num]][i2] <- 0
    }
    if (nb.model > nb.model0)
      res[(nb.model0 + 1):nb.model, c("ll", "df", "criterion")] <- NA
    res <- cbind(G == 1, res)
    if (nb.model == 1)
      res[, c("ll", "df", "criterion")] <- NULL
    row.names(res) <- 1:nrow(res)
    return(list(model = lm.min, res = res, cov0 = cov0))
  }


# additional function -----------------------------------------------------


modelFromSelection <- function(Y,X,selection=NULL){
  if(!is.data.frame(Y)){
    Y <- as.data.frame(Y)
  }
  if(!is.data.frame(X)){
    X <- as.data.frame(X)
  }

  m = ncol(Y)
  lm.list=list()
  for(i in 1:m){
    data = cbind(y = matrix(Y[,i],ncol=1),X)
    formula = paste0("y ~ 1 ")
    if(m==1){list.c = names(selection)[selection==1]}else{list.c = colnames(selection)[selection[i,]==1]}
    if(length(list.c)>0){
      formula = paste0(formula," + ",paste0(list.c,collapse=" + "))
    }
    lm.sel <- NULL
    eval(parse(text=paste0("lm.sel <- lm(",formula,",data=data)")))
    lm.list <- append(lm.list,list(lm.sel))
  }

  if(m==1){
    return(lm.sel)
  }else{
    return(lm.list)
  }
}

## Cette fonction prends en entrée une matricede covariables et d'individus X, Y tel que :
##     -  Y est un vecteur donc pour UN paramètre avec en ligne les différents individus
##     - Y est une matrice avec sur chaque colonne le paramètre
##
## Sigma est la matrice de covariance des paramètres considérés
## X est la matrice des covariables
##
## A ce niveau là, rien n'est BLANCHI ni STANDARDIZE, X et Y sont des matrices quantitatives.
## cov0 correspond au colonne de X a exclure de la sélection (doit être du même type)

updateCov0 <- function(Y, eta.list, x, p.max=1, covFix = NULL,
                       pen.coef = NULL, pw.list = NULL, cov0.list = NULL){
  param = names(cov0.list)
  n.param=length(param)
  t.param = setdiff(colnames(Y),c("id","rep"))
  for(j in 1:n.param){
    cov0.list[[param[j]]] <- upd(param[j],Y[,t.param[j],drop=FALSE],eta.list[[param[j]]],x,p.max,covFix,pen.coef,pw.list[[param[j]]],cov0.list[[param[j]]])
  }
  return(cov0.list)
}

upd <- function(ny,y,eta,x,p.max,covFix,pen.coef,pw,cov0){
  N <- length(unique(x$id))
  nrep <- nrow(y)/N
  if (p.max <= 1) {
    nx <- setdiff(names(x), c("id", "rep"))
    yc <- rowMeans(matrix(y[[1]], nrow = N))
    xc <- x
    lm0 <- lm(yc ~ 1)
    nxc <- setdiff(nx, cov0)
    pjc <- NULL
    for (nc in nxc) {
      lmc <- lm(yc ~ xc[[nc]])
      pc <- signif(anova(lm0, lmc)$`Pr(>F)`[2], 4)
      pjc <- c(pjc, pc)
    }
    pjc <- Rsmlx.p.weight(pjc, pw[nxc], pen.coef)
    names(pjc) <- nxc
    if (!is.null(eta)) {
      etac <- rowMeans(matrix(eta[[1]], nrow = N))
      lm0 <- lm(etac ~ 1)
      pjec <- NULL
      for (nc in nxc) {
        lmc <- lm(etac ~ xc[[nc]])
        pc <- signif(anova(lm0, lmc)$`Pr(>F)`[2], 4)
        pjec <- c(pjec, pc)
      }
      pjec <- Rsmlx.p.weight(pjec, pw[nxc], pen.coef)
      names(pjec) <- nxc
      pjc <- pmin(pjc, pjec, na.rm = TRUE)

    }
    pjc[names(which(mlx.getIndividualParameterModel()$covariateModel[[ny]]))] <- 0
    list.c <- which(pjc > p.max)
    cov0 <- c(cov0, nxc[list.c])
  }
  return(cov0)
}

