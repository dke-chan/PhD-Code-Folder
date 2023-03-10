library(ascr)
source("./01-aux-funcs.R")
sim.theta = readRDS(file = "./model_fit.RDS")
survey.design = readRDS(file = "./surveyDesigns.RDS")
seeds = read.csv(file = "./seeds.csv")

which.design <- 1

run.id <- 2
cluster.work = subset(seeds, run.ids == run.id)
toFit = expand.grid(detfn = c("ss"), aux = c(FALSE), stringsAsFactors = FALSE)

for (i in 1:nrow(cluster.work)) {
  set.seed(cluster.work$seed[i])

  ## Sim the data
  y <- try(sim.capt.wrap("fascr", survey.design[[which.design]]))

  while (inherits(y, "try-error")) {
    y <- try(sim.capt.wrap("fascr", survey.design[[which.design]]))
  }

  cat(cluster.work$seed[i], "\n")

  for (j in 1:nrow(toFit)) {
    cat("    ",toFit$detfn[j], toFit$aux[j], "\n")

    ## Setup the capture history matrix properly
    capt <- y[[1]]
    if (toFit$detfn[j] != "ss") {
      capt$ss <- NULL
    }
    if (!toFit$aux[j]) {
      capt$toa <- NULL
    }

    ##
    star.fit <- try(fit.ascr(capt, survey.design[[which.design]], y[[2]], detfn = toFit$detfn[j],
                             ss.opts = sim.theta[["fascr"]]$ss.opts,
                             fix = list(b0.ss = 150)))
                             # sv = list(D = sim.theta[["fascr"]]$theta$D)))
    if (!inherits(star.fit, what = "try-error")) {
      cat("    logLik SV:", logLik(star.fit), "\n")
    } else {
      cat("    logLik SV: failed \n")
    }

    star.fit.2 <- try(fit.ascr(capt, survey.design[[which.design]], y[[2]], detfn = toFit$detfn[j],
                             ss.opts = sim.theta[["fascr"]]$ss.opts,
                             fix = list(b0.ss = 175)))

    if (!inherits(star.fit.2, what = "try-error")) {
      cat("    logLik no SV:", logLik(star.fit.2), "\n")
    } else {
      cat("    logLik SV: failed \n")
    }
  }


}


