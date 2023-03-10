# Setup data etc. ---------------------------------------------------------

load("./whales.RData") # load("./data/whales.RData")
source("./00-aux-funcs.R") # source("./code/whales/00-aux-funcs.R")
hpp.fits <- read.csv("./hpp-fits.csv", row.names = 1) # hpp.fits <- read.csv("./code/whales/hpp-fits.csv", row.names = 1)

library(ascr); library(mgcv)

cali.Capt <- list(bincapt = as.matrix(cali.CaptHist), toa = as.matrix(cali.ToA))
cali.covariates <- data.frame(depth = cali.Mask.z[cali.Mask.z < 0], coast = cali.Mask.minCoastDist)

# Load in the model the core will fit -------------------------------------

run.id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")) # 146
model.input <- read.csv("./scr-to-fit.csv")[run.id, ] # model.input <- read.csv("./code/whales/scr-to-fit.csv")[run.id, ]

# Fit the model without estimating the Hessian ----------------------------

mod.fit.1 <- coreWorkerFunc(model.input)
mod.fit.1$value$args$hess <- TRUE

# Fit the model with starting values (if necessary) and estimate the Hessian ------

mod.fit.2 <- tryCatch(
        withWarnings(do.call("fit.ascr", mod.fit.1$value$args)),
        error = function (e) e)

# Save results.. ----------------------------------------------------------

toSave <- list()
rect.covariate.outer <- try(mod.fit.2$value$scale.covs(covariates.traps.df)) # Naming this variable like this because of rectDensityEstimate()..

toSave$AIC <- try(AIC(mod.fit.2$value))
toSave$numOfWarn <- length(mod.fit.2$warnings)
toSave$coef <- try(get.par(mod.fit.2$value))
toSave$vcov <- try(vcov(mod.fit.2$value))
toSave$rectDesignMatrix <- try(predict(gam(G = mod.fit.2$value$fgam), newdata = rect.covariate.outer, type = "lpmatrix"))
toSave$warnMsg <- mod.fit.2$warnings
toSave$maxgrad <- try(mod.fit.2$value$maxgrad)
toSave$rectDensity <- try(rectDensityEstimate(mod.fit.2$value))

saveRDS(toSave, paste0("./full-grid-results/", run.id, ".rds")) # saveRDS(toSave, paste0("./code/whales/full-grid-results/", run.id, ".rds"))
