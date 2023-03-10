# Setup data etc. ---------------------------------------------------------

load("./data/whales.RData")
source("./code/whales/00-aux-funcs.R")

all.fits <- list.files("./code/whales/full-grid-results/")

all.fits.id <- all.fits |>
  gsub(pattern = "([0-9]*)\\.rds", replacement = "\\1") |>
  as.numeric()

all.fits.ref <- data.frame(
  filePath = paste0("./code/whales/full-grid-results/", all.fits),
  orderID = all.fits.id
)

all.fits.ref <- all.fits.ref[order(all.fits.ref$orderID), ]

library(ascr); library(mgcv)

results.Master <- vector("list", nrow(all.fits.ref))

for (i in 1:length(results.Master)) {
  results.Master[[i]] <- readRDS(all.fits.ref$filePath[i])
}

# -------------------------------------------------------------------------

full.set <- read.csv("./code/whales/scr-to-fit.csv")

full.set$modelNum <- 1:nrow(full.set)

full.set$BIC <- lapply(results.Master, \(x) {
  toReturn <- try(x$AIC)
  if (inherits(toReturn, "try-error")) {
    return (0)
  }

  ll <- (toReturn - 2 * length(x$coef - 1)) / -2

  return ((length(x$coef) - 1) * log(nrow(cali.CaptHist)) - 2 * ll)}) |>
  unlist()

full.set$n.par <- lapply(results.Master, \(x) {
  length(x$coef) - 1 }) |>
  unlist()

full.set$AIC <- lapply(results.Master, \(x) {
  toReturn <- try(x$AIC)
  if (inherits(toReturn, "try-error")) {
    return (0)
  }
  return (toReturn)}) |>
  unlist()

full.set$noWarn <- lapply(results.Master, \(x) {try(x$numOfWarn)}) |>
  unlist()

full.set$maxgrad <- 0
full.set$maxgrad[which(full.set$AIC != 0)] <- lapply(results.Master, \(x) {try(x$maxgrad)}) |>
  unlist()

full.set$rectEst <- lapply(results.Master, \(x) {
  toReturn <- try(x$rectDensity[["estimate"]])
  if (inherits(toReturn, "try-error")) {
    return (0)
  }
  return (toReturn)}) |>
  unlist()

full.set$rectStdErr <- lapply(results.Master, \(x) {
  toReturn <- try(x$rectDensity[["standard.error"]])
  if (inherits(toReturn, "try-error")) {
    return (0)
  }
  return (toReturn)}) |>
  unlist()

full.set$rectLwr <- lapply(results.Master, \(x) {
  toReturn <- try(x$rectDensity[["wald.interval"]][1])
  if (inherits(toReturn, "try-error")) {
    return (0)
  }
  return (toReturn)}) |>
  unlist()

full.set$rectUpr <- lapply(results.Master, \(x) {
  toReturn <- try(x$rectDensity[["wald.interval"]][2])
  if (inherits(toReturn, "try-error")) {
    return (0)
  }
  return (toReturn)}) |>
  unlist()

full.set$det <- lapply(results.Master, \(x) {
  toReturn <- try(det(x$vcov))
  if (inherits(toReturn, "try-error")) {
    return (NA)
  }
  return (toReturn)})

# What were the best three models? ----------------------------------------

numOfModels <- 3
results.sorted <- subset(full.set, noWarn == 0 & AIC != 0 & !is.nan(rectStdErr)) |>
  {\(x){
    ix <- sort(x$BIC, index.return = TRUE)$ix
    x[ix, ]
  }}()

toTabulate <- results.sorted[1:numOfModels, ]

missingDetFn <- c("hn", "hn.f", "hr", "th", "hhn") %in% toTabulate$detfn |>
  {\(x) {c("hn", "hn.f", "hr", "th", "hhn")[!x]}}()

for (detfn in missingDetFn) {
  toTabulate <- rbind(
    toTabulate,
    results.sorted[which(results.sorted$detfn == detfn)[1], ]
  )
}

saveRDS(toTabulate, "./code/whales/ihpp-fitted.rds")
saveRDS(full.set, "./code/whales/ihpp-all.rds")

# -------------------------------------------------------------------------

library(parallel)
cl.cores <- makeCluster(4)

cali.Capt <- list(bincapt = as.matrix(cali.CaptHist), toa = as.matrix(cali.ToA))
cali.covariates <- data.frame(depth = cali.Mask.z[cali.Mask.z < 0], coast = cali.Mask.minCoastDist)

clusterEvalQ(cl.cores, {
  library(ascr)
  load("./data/whales.RData")
  hpp.fits <- read.csv("./code/whales/hpp-fits.csv", row.names = 1)
})
clusterExport(cl.cores, varlist = c(
  "model_fit", "withWarnings", "cali.Capt", "cali.covariates", "full.set", "coreWorkerFunc"
))

toTabulate.models <- parLapplyLB(cl.cores, toTabulate$modelNum, function (modelIndex) {
  x <- coreWorkerFunc(full.set[modelIndex, ])
  x$value$args$hess <- TRUE

  return (try(do.call("fit.ascr", x$value$args)))})

stopCluster(cl.cores)
remove(cl.cores)

for (i in 1:nrow(toTabulate)) {
  saveRDS(toTabulate.models[[i]], paste0("./code/whales/ihpp-model-", toTabulate$modelNum[i], ".rds"))
}
