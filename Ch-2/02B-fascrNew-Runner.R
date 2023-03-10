library(ascr)
source("./01-aux-funcs.R")
sim.theta = readRDS(file = "./model_fit.RDS")
survey.design = readRDS(file = "./surveyDesigns.RDS")
seeds = read.csv(file = "./seeds.csv")

# -------------------------------------------------------------------------

run.id <- 2 # as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
which.design <- 1

# -------------------------------------------------------------------------

cluster.work = subset(seeds, run.ids == run.id)
toFit = expand.grid(detfn = c("hn", "hhn", "hr", "ss"), aux = c(TRUE, FALSE), stringsAsFactors = FALSE)

# -------------------------------------------------------------------------

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

    alreadyDone = read.csv(paste0("./fascrNew_", which.design, ".csv")) |>
      subset(seed == cluster.work$seed[i] & detfn == toFit$detfn[j] & aux == toFit$aux[j]) |>
      nrow()

    if (alreadyDone > 0) {
      next
    }

    ## Setup the capture history matrix properly
    capt <- y[[1]]
    if (toFit$detfn[j] != "ss") {
      capt$ss <- NULL
    }
    if (!toFit$aux[j]) {
      capt$toa <- NULL
    }

    fix.list = NULL
    ## Fix g0=1 for 2-by-1 when using a "hn" function without aux data
    if (toFit$detfn[j] == "hn" & !toFit$aux[j] & which.design == 1) {
      fix.list = list(g0 = 1)
    }

    ## Fix g0=1 for 2-by-1, 3-by-1, and triangle when using a "hr" function without aux data
      ## You may need to
    if (toFit$detfn[j] == "hr" & !toFit$aux[j] & (which.design == 1 | which.design == 2 | which.design == 4)) {
      fix.list = list(g0 = 1)
    }

    ##
    star.fit <- try(fit.ascr(capt, survey.design[[which.design]], y[[2]], detfn = toFit$detfn[j],
                             ss.opts = sim.theta[["fascr"]]$ss.opts,
                             sv = list(D = sim.theta[["fascr"]]$theta$D), fix = fix.list))

    ##
    if (inherits(star.fit, "ascr")) {
      theta <- c(get.par(star.fit, pars = "fitted"), star.fit$loglik)
      names(theta)[length(theta)] <- "loglik"

      toWrite <- data.frame(theta = theta, theta.name = names(theta), survey = "fascr",
                            detfn = toFit$detfn[j], seed = cluster.work$seed[i], aux = toFit$aux[j])
    } else {
      toWrite <- data.frame(theta = NA, theta.name = NA, survey = "fascr",
                            detfn = toFit$detfn[j], seed = cluster.work$seed[i], aux = toFit$aux[j])
    }

    write.table(toWrite, file = paste0("./fascrNew_", which.design, ".csv"), # "./code/misspecified/fascr_72.csv",
                append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
  }


}


