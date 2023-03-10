library(ascr); library(pracma); library(tibble); library(CircStats); library(parallel)
source("./code/sonobuoys/00-aux-funcs.R")
load("./data/sonobuoys.RData")
logistic.like.models <- readRDS("./code/sonobuoys/logistic-like-models.RDS")

# von Mises Truth ---------------------------------------------------------
act.bearing <- seaTrial.df$act.bearing.rad
obs.bearing <- seaTrial.df$obs.bearing.rad
which.Indices <- seaTrial.df$Detected != 0

bearing.err <- (act.bearing - obs.bearing)[which.Indices]

kappa.theta <- CircStats::est.kappa(bearing.err)

# -------------------------------------------------------------------------

numOfSims <- 1000
detfn.form <- c("hn", "hhn", "hr", "th")
detfn.theta <- lapply(detfn.form, \(x) {
  logistic.like.models$est[logistic.like.models$detfn == x]
}) |> setNames(detfn.form)

set.seed(2787484)
seeds <- sample(1:1e7, numOfSims * length(detfn.form), replace = FALSE)

simMaster <- tibble(detfn = rep(detfn.form, each = numOfSims), theta = rep(detfn.theta, each = numOfSims), kappa = kappa.theta, seed = seeds)

fileNames <- paste0("./code/sonobuoys/", detfn.form, "Sim.csv")

if (any(!file.exists(fileNames))) {
  for (i in fileNames) {write.table(data.frame(output = numeric(0), whichOutput = character(0), model = character(0), detfn = character(0), seed = numeric(0), seed.index = numeric(0)), file = i, row.names = FALSE, sep =",")}

  cluster <- makeCluster(7)

  clusterEvalQ(cluster, {
    library(ascr); library(pracma); library(tibble); library(CircStats); library(parallel)
    source("./code/sonobuoys/00-aux-funcs.R")
    load("./data/sonobuoys.RData")
  })

  clusterExport(cluster, c("simMaster"))

  hide <- clusterApplyLB(cluster, 1:nrow(simMaster), \(x) {
    set.seed(simMaster$seed[x])

    ascr.fit <- NULL
    mrds.fit <- NULL
    logi.fit <- NULL
    failed.it <- -1
    detfn.param <- switch(simMaster$detfn[x], hn = c("g0", "sigma"), hr = c("g0", "sigma", "z"),
                          th = c("shape", "scale"), hhn = c("lambda0", "sigma"))
    detfn.f <- switch(simMaster$detfn[x], hn = halfnormal, hr = hazardRate,
                      th = threshold, hhn = hazardHalfnormal)
    start.val <- as.list(simMaster$theta[[x]]) |>
      setNames(detfn.param)

    while (!(inherits(ascr.fit, "ascr") & inherits(mrds.fit, "ascr") & inherits(logi.fit, "list"))) {
      failed.it <- failed.it + 1
      data.star <- simNewCaptures(simMaster$detfn[x], simMaster$theta[[x]], simMaster$kappa[x])

      ascr.fit <- tryCatch(
        fit.ascr(capt = lapply(data.star$capt.obj, \(x) {x$mrds <- NULL; return(x)}), traps = data.star$traps.obj, mask = data.star$mask.obj, detfn = simMaster$detfn[x], sv = start.val),
        error = function (cnd) {NA},
        warning = function (cnd) {NA}
      )

      if (inherits(ascr.fit, "ascr")) {
        mrds.fit <- tryCatch(
          fit.ascr(capt = data.star$capt.obj, traps = data.star$traps.obj, mask = data.star$mask.obj, detfn = simMaster$detfn[x], sv = start.val),
          error = function (cnd) {NA},
          warning = function (cnd) {NA}
        )

        if (inherits(mrds.fit, "ascr")) {
          bernoulli.trial <- as.vector(data.star$capt.obj.full)
          bernoulli.distance <- as.vector(full.distance.mat * 1000)

          logi.fit <- tryCatch(
            optim(c(1, log(1000), log(10))[1:length(start.val)], bernoulliLL, probFunc = detfn.f, y = bernoulli.trial, x = bernoulli.distance, hessian = TRUE, method = "L-BFGS-B"),
            error = function (cnd) {NA},
            warning = function (cnd) {NA}
          )
        }
      }
    }

    ascr.toReturn <- modelResPrinter(ascr.fit, "ascr")[, c("est", "est.name", "spec", "detfn")]
    mrds.toReturn <- modelResPrinter(mrds.fit, "mrds")[, c("est", "est.name", "spec", "detfn")]

    to.exp <- which(detfn.param %in% c("sigma", "lambda0", "scale", "z"))
    to.plogis <- which(detfn.param %in% c("g0"))
    logi.toReturn <- data.frame(est = logi.fit$par, est.name = detfn.param, spec = "logi", detfn = simMaster$detfn[x])
    logi.toReturn$est[to.exp] <- exp(logi.fit$par[to.exp])
    logi.toReturn$est[to.plogis] <- plogis(logi.fit$par[to.plogis])

    L2Norms <- c(euclideanNorm(ascr.toReturn$est[ascr.toReturn$est.name %in% detfn.param], logi.toReturn$est, detfn.f, link = FALSE),
                 euclideanNorm(mrds.toReturn$est[mrds.toReturn$est.name %in% detfn.param], logi.toReturn$est, detfn.f, link = FALSE),
                 euclideanNorm(ascr.toReturn$est[ascr.toReturn$est.name %in% detfn.param], mrds.toReturn$est[mrds.toReturn$est.name %in% detfn.param], detfn.f, link = FALSE))


    toPrint <- rbind(ascr.toReturn, mrds.toReturn, logi.toReturn,
                     data.frame(est = c(L2Norms, failed.it, length(data.star$capt.obj)), est.name = c("A-L", "M-L", "A-M", "failed.it", "n.detect"),
                                spec = "stats", detfn = simMaster$detfn[x]))

    toPrint$seed <- simMaster$seed[x]
    toPrint$seed.index <- x

    colnames(toPrint)[1:3] <- c("output", "whichOutput", "model")

    write.table(x = toPrint, file = paste0("./code/sonobuoys/", simMaster$detfn[x], "Sim.csv"), append = TRUE,
                sep = ",", row.names = FALSE, col.names = FALSE)

    return (x)
  })

  stopCluster(cluster)
}
