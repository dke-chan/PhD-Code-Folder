model_fit <- function(detfn, covariate, spline.basis, spline.k, ...) {
  if (covariate == "none") {
    toFormula <- "~ 1"
  } else if (covariate != "both") {
    if (spline.k[1] == 1) {
      toFormula <- paste0("~", covariate)
    } else if (spline.k[1] == 2) {
      toFormula <- paste0("~", covariate, " + I(", covariate, "^2)")
    } else {
      toFormula <- paste0("~s(", covariate, ", k = ", spline.k[1],",bs = \"", spline.basis[1],"\")")
    }
  } else {
    if (spline.k[1] == 1) {
      toFormula.A <- paste0("~ coast")
    } else if (spline.k[1] == 2) {
      toFormula.A <- paste0("~ coast + I(coast^2)")
    } else {
      toFormula.A <- paste0("~s(coast, k = ", spline.k[1],",bs = \"", spline.basis[1],"\")")
    }

    if (spline.k[2] == 1) {
      toFormula.B <- paste0("+ depth")
    } else if (spline.k[2] == 2) {
      toFormula.B <- paste0("+ depth + I(depth^2)")
    } else {
      toFormula.B <-  paste0("+ s(depth, k = ", spline.k[2],", bs = \"", spline.basis[1],"\")")
    }

    toFormula <- paste0(toFormula.A, toFormula.B)
  }

  fix.list <- list()
  if (detfn == "hn.f") {
    detfn <- "hn"
    fix.list <- list(g0 = 1)
  }

  # print(toFormula)

  fit.ascr(capt = cali.Capt,
           traps = cali.Traps,
           mask = cali.Mask.Sea,
           detfn = detfn,
           sound.speed = 1500,
           ihd.opts = list(model = as.formula(toFormula),
                           covariates = cali.covariates),
           fix = fix.list,
           ...)
}

# -------------------------------------------------------------------------

withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}

# -------------------------------------------------------------------------

coreWorkerFunc <- function(input) {

  model.built <- tryCatch(
    withWarnings(model_fit(input$detfn, input$covariate.input, input$basis, as.numeric(input[, c("k.1", "k.2")]), hess = FALSE)),
    error = function (e) e)

  if (inherits(model.built$value, "ascr") & length(model.built$warnings) > 0) {
    ## Attempt to converge at the fitted values
    warn.built <- tryCatch(
      withWarnings(model_fit(input$detfn, input$covariate.input, input$basis, as.numeric(input[, c("k.1", "k.2")]), hess = FALSE,
                             sv = get.par(model.built$value, "fitted", as.list = TRUE))),
      error = function (e) e)

    if (inherits(warn.built$value, "ascr") & length(warn.built$warnings) > 0) {
      get.base.par <- get.par(model.built$value, "fitted", as.list = TRUE)
      get.base.par[grepl("D", names(get.base.par))] <- NULL
      get.base.par$sigma.toa <- subset(hpp.fits, detfn == input$detfn)$sigma.toa
      get.base.par[!grepl("sigma.toa", names(get.base.par))] <- subset(hpp.fits, detfn == input$detfn)[paste0("theta.", 1:3)] |>
        {\(x){x[!is.na(x)]}}()

      warn.built <- tryCatch(
        withWarnings(model_fit(input$detfn, input$covariate.input, input$basis, as.numeric(input[, c("k.1", "k.2")]), hess = FALSE,
                               sv = get.base.par)),
        error = function (e) e)
    }

    if (!inherits(warn.built, "simpleError")) {
      return (warn.built)
    }
  } else if (inherits(model.built, "simpleError")) {
    ##
    error.built <- tryCatch(
      withWarnings(model_fit(input$detfn, input$covariate.input, input$basis, as.numeric(input[, c("k.1", "k.2")]), hess = FALSE,
                             sv = get.par(model.built$value, "fitted", as.list = TRUE))),
      error = function (e) e)

    if (inherits(error.built, "simpleError")) {
      return (error.built)
    }
  }

  return (model.built)
}

# -------------------------------------------------------------------------

simCapt <- function(ascrObj, detDistance = 15000) {
  density.Surf <- predict(ascrObj)
  traps <- get.traps(ascrObj)
  mask.points <- get.mask(ascrObj)
  ascrObj.par <- get.par(ascrObj, "fitted", as.list = TRUE)

  #
  xlim <- range(traps[, 1]) + c(-1, 1) * detDistance
  ylim <- range(traps[, 2]) + c(-1, 1) * detDistance

  xbound <- xlim[1] <= mask.points[, 1] & mask.points[, 1] <= xlim[2]
  ybound <- ylim[1] <= mask.points[, 2] & mask.points[, 2] <= ylim[2]

  mask.points.cut <- convert.mask(mask.points[xbound & ybound, ])

  #
  call.star <- secr::sim.popn(D = density.Surf[[1]][xbound & ybound], core = mask.points.cut, buffer = 0, model2D = "IHP") |>
    as.matrix()

  #
  n.calls <- nrow(call.star)
  dists <- distances(call.star, traps)
  n.traps <- nrow(traps)

  det.probs <- ascr:::calc.detfn(dists, ascrObj$args$detfn, ascrObj.par[ascrObj$detpars])

  full.bin.capt <- matrix(as.numeric(runif(n.calls * n.traps) < det.probs), nrow = n.calls, ncol = n.traps)
  captures <- which(apply(full.bin.capt, 1, sum) > 0)
  bin.capt <- full.bin.capt[captures, , drop = FALSE]

  n.dets <- sum(bin.capt)
  capt.popn <- call.star[captures, ]
  capt.dists <- dists[captures, ]

  toa.capt <- capt.dists/ascrObj$args$sound.speed * bin.capt
  toa.capt[bin.capt == 1] <- toa.capt[bin.capt == 1] + rnorm(n.dets, sd = ascrObj.par$sigma.toa)

  return(list(bincapt = bin.capt, toa = toa.capt))
}

# -------------------------------------------------------------------------

boostrapFit <- function(whichSeed, whichSource, whichThresh = 16000) {
  toReturn <- list(numTries = 0)
  set.seed(whichSeed)

  args <- whichSource$args
  args$capt <- simCapt(whichSource, whichThresh)
  args$hess <- FALSE

  ascr.boot <- tryCatch(withWarnings(do.call("fit.ascr", args)),
                        error = function (e) e)

  while (inherits(ascr.boot$value, "ascr") & length(ascr.boot$warnings) > 0) {
    toReturn$numTries <- toReturn$numTries + 1
    args$capt <- simCapt(whichSource, whichThresh)
    ascr.boot <- tryCatch(withWarnings(do.call("fit.ascr", args)),
                          error = function (e) e)
  }

  toReturn$coef <- get.par(ascr.boot$value)
  toReturn$warnings <- ascr.boot$warnings
  toReturn$rect <- predict(ascr.boot$value, newdata = covariates.traps.df)

  return (toReturn)
}

# -------------------------------------------------------------------------

rectDensityEstimate <- function(whichModel, style.two = TRUE) {
  thetaHat <- coef(whichModel)
  SigmaHat <- vcov(whichModel)
  rect.covariate <- whichModel$scale.covs(covariates.traps.df)
  designMatrix <- predict(gam(G = whichModel$fgam), newdata = rect.covariate, type = "lpmatrix")

  betaIndices <- grep(pattern = "D.", x = names(thetaHat))
  scrIndices <- grep(pattern = "D.", x = names(thetaHat), invert = TRUE)

  if (style.two) {
    exp_summation <- designMatrix %*% thetaHat[betaIndices] |> exp()
    numerator <- colSums(designMatrix * matrix(rep(exp_summation, times = ncol(designMatrix)), nrow = nrow(designMatrix), ncol = ncol(designMatrix)))

    jacob <- c(numerator / sum(exp_summation), rep(0, length(scrIndices)))
    std.err <- t(jacob) %*% SigmaHat %*% jacob |> as.numeric() |> sqrt()

    log_est <- predict(whichModel, newdata = covariates.traps.df) |> mean() |> log()

    return (list(
      estimate = exp(log_est),
      standard.error = std.err,
      wald.interval = exp(log_est + c(-1, 1) * qnorm(0.975) * std.err)
    ))
  } else {
    jacob <- c(colSums(designMatrix) / nrow(rect.covariate), rep(0, length(scrIndices)))
    std.err <- t(jacob) %*% SigmaHat %*% jacob |> as.numeric() |> sqrt()

    est <- predict(whichModel, newdata = covariates.traps.df, use.log = TRUE) |> mean()
    std.err <- matrix(jacob, nrow = 1) %*% SigmaHat %*% matrix(jacob, ncol = 1) |> as.numeric() |> sqrt()

    return (list(
      estimate = exp(est),
      standard.error = std.err,
      wald.interval = exp(est + c(-1, 1) * qnorm(0.975) * std.err)
    ))
  }
}
