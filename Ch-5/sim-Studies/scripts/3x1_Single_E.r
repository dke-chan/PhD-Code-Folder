#
library(moveLLR)
library(spatstat.geom, include.only = c("pairdist", "crossdist"))
library(spatstat.utils, include.only = "matrowany")
library(CircStats, include.only = "rvm")
library(mvtnorm, include.only = "rmvnorm")
library(dplyr)
library(ascr)
library(secr)
source("./simul.r")

#
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#
arrayRef <- read.csv("./arrayRef.csv") |>
  subset(arrayID == id)
design <- "3x1_Single_E"

#
listeningPosts <- expand.grid(x = c(-1000, 0, 1000), y = c(0))
nOccasion <- 3
sigmaMove <- 1000

#
std.start <- list(D = 0.5 / 100 * 2/3 * nOccasion, lambda0 = 5, sigma = 750, kappa = 10)
stack.start <- c(log(0.5 / 100), log(5), log(750), qlogis(2/3), log(10))
disperse.start <- c(log(0.5 / 100), log(5), log(750), qlogis(2/3), log(sigmaMove), log(10))

#
aBounds <- data.frame(x = range(listeningPosts$x) + 15000 * c(-1, 1), y = range(listeningPosts$y) + 15000 * c(-1, 1))
mask <- ascr::create.mask(as.matrix(listeningPosts), buffer = 9000)
maskDists <- ascr::distances(mask, as.matrix(listeningPosts))
expectedBearings <- t(ascr::bearings(as.matrix(listeningPosts), mask))

#
for (run in arrayRef$seed) {
  set.seed(run)

  full.df <- simul_data(
    density = 0.5, survey.region = aBounds, listening.post = listeningPosts,
    detfn.theta = c(lambda0 = 5, sigma = 750), kappa = 10,
    movement.theta = c(sigma_11 = sigmaMove, sigma_12 = 0, sigma_21 = 0, sigma_22 = sigmaMove)^2,
    call.Process = "Binomial", call.Process.Theta = c(nOccasion, 2/3)
  )

  # Create the capture history and their bearings from `simul.frame`
  bin.capt <- full.df %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    select(starts_with("detect.")) %>%
    as.matrix()

  bearings <- full.df %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    select(starts_with("ob.bearing")) %>%
    as.matrix()

  bearings[, 1:nrow(listeningPosts)] <- bearings[, 1:nrow(listeningPosts)] * bin.capt[, 1:nrow(listeningPosts)]

  # Pull the group identifiers out
  animal.id <- full.df %>%
    mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
    filter(detected > 0) %>%
    pull(groupID) %>%
    factor() %>%
    as.numeric() %>%
    {. - 1}

  #
  cat("Seed:", run, "\n")

  # Fit the models
  std.fit <- try(fit.ascr(
    capt = list(bincapt = bin.capt, bearing = bearings), traps = as.matrix(listeningPosts), mask = mask,
    detfn = "hhn", sv = std.start, trace = TRUE
  ))

  cat("---------------\n")

  stack.fit <- try(optim(
    par = stack.start, fn = fit_ascrStack, bin_capt = bin.capt, traps = as.matrix(listeningPosts), mask = mask,
    observed_bearings = bearings, expected_bearing = expectedBearings, mask_dists = maskDists, pixel_area = attr(mask, "area"),
    ID = animal.id, survey_length = nOccasion, hessian = TRUE, method = "L-BFGS-B", bernoulli = TRUE
  ))

  cat("---------------\n")

  disperse.fit <- try(optim(
    par = disperse.start, fn = fit_ascrDisperse, bin_capt = bin.capt, traps = as.matrix(listeningPosts), mask = mask,
    observed_bearings = bearings, expected_bearing = expectedBearings, mask_dists = maskDists, pixel_area = attr(mask, "area"),
    ID = animal.id, survey_length = nOccasion, hessian = TRUE, method = "Nelder-Mead", trace = TRUE
  ))

  #
  cat("===============\n")

  # Save results
  if (!inherits(std.fit, "try-error")) {
    results <- data.frame(seed = run, theta = get.par(std.fit, pars = "fitted"), se = sqrt(diag(vcov(std.fit, pars = "linked"))), model = "ascr", design = design)
  } else {
    results <- data.frame(seed = run, theta = rep(NA, length(std.start)), se = rep(NA, length(std.start)), model = "ascr", design = design)
  }
  write.table(results, file = paste0("./results/", design, ".csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

  if (!inherits(stack.fit, "try-error")) {
    se.toReturn <- try(sqrt(diag(solve(stack.fit$hessian))))
    if (inherits(se.toReturn, "try-error")) {
      se.toReturn <- c(rep(NA, length(stack.start)), NA)
    } else {
      se.toReturn <- c(se.toReturn, solve(stack.fit$hessian)[1, 4])
    }
    results <- data.frame(seed = run, theta = c(stack.fit$par, NA), se = se.toReturn, model = "ascrStack", design = design)
  } else {
    results <- data.frame(seed = run, theta = c(rep(NA, length(stack.start)), NA), se = rep(NA, length(stack.start) + 1), model = "ascrStack", design = design)
  }
  write.table(results, file = paste0("./results/", design, ".csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")

  if (!inherits(disperse.fit, "try-error")) {
    se.toReturn <- try(sqrt(diag(solve(disperse.fit$hessian))))
    if (inherits(se.toReturn, "try-error")) {
      se.toReturn <- c(rep(NA, length(disperse.start)), NA)
    } else {
      se.toReturn <- c(se.toReturn, solve(disperse.fit$hessian)[1, 4])
    }
    results <- data.frame(seed = run, theta = c(disperse.fit$par, NA), se = se.toReturn, model = "ascrDisperse", design = design)
  } else {
    results <- data.frame(seed = run, theta = c(rep(NA, length(disperse.start)), NA), se = rep(NA, length(disperse.start) + 1), model = "ascrDisperse", design = design)
  }
  write.table(results, file = paste0("./results/", design, ".csv"), append = TRUE, row.names = FALSE, col.names = FALSE, sep = ",")
}
