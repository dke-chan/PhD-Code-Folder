load("./data/whales.RData")
source("./code/whales/00-aux-funcs.R")
library(sp); library(spatstat); library(mgcv)

# -------------------------------------------------------------------------

scrDateRange <- read.table("Data/capt-toa.txt") |>
  apply(2, \(x) {as.POSIXct(ifelse(x > 0, (x - 719529) * 86400, 0), origin = "1970-01-01")}) |>
  as.numeric() |>
  {\(x){x[x > 0]}}() |>
  range() |>
  as.POSIXct(origin = "1970-01-01")

# Sightings
sightings.full <- read.csv("Data/obis_seamap_custom_5d33f23f136fa_20190721_011256_points_csv.csv", stringsAsFactors = FALSE)
sightings.full$date_time_2 <- as.POSIXct(sightings.full$date_time, format = "%Y-%m-%d %H:%M:%S")
sightings.df <- subset(sightings.full, scrDateRange[1] <= date_time_2 & date_time_2 <= scrDateRange[2])

# Camera location
# 36 26' 23.61" N 121 55' 20.5" W from http://seamap.env.duke.edu/dataset/1862
camera.loc <- data.frame(latitude = 36.439892, longitude = -121.922361)

# Coordinates system change
coordinates(sightings.full) <- c("longitude", "latitude")
proj4string(sightings.full) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
sightings.full <- spTransform(sightings.full, CRS("+proj=utm +zone=10 +datum=WGS84"))

coordinates(sightings.df) <- c("longitude", "latitude")
proj4string(sightings.df) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
sightings.df <- spTransform(sightings.df, CRS("+proj=utm +zone=10 +datum=WGS84"))

coordinates(camera.loc) <- c("longitude", "latitude")
proj4string(camera.loc) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
camera.loc <- spTransform(camera.loc, CRS("+proj=utm +zone=10 +datum=WGS84"))

# Making mask.owin
# whichPoints <- chull(sightings.full@coords[, 1], sightings.full@coords[, 2])
# mask.owin <- owin(poly = sightings.full@coords[rev(whichPoints), ])

# Converting the sightings to a spatstat object
sightings.pp <- ppp(sightings.df@coords[, 1], sightings.df@coords[, 2],
                    xrange = range(sightings.df@coords[, 1]), yrange = range(sightings.df@coords[, 2]))
                    # mask.owin)

# Spatial covariates
sightings.covariates <- cbind(cali.Mask.Sea, sapply(cali.covariates, scale)) |>
  as.data.frame() |>
  as.im()

# -------------------------------------------------------------------------

ppm.full.set <- expand.grid(k.1 = 0:4,
                            k.2 = 0:4,
                            whichInput = c("coast", "depth", "both"),
                            stringsAsFactors = FALSE) |>
  subset(!(whichInput == "both" & (k.1 == 0 | k.2 == 0))) |>
  subset(!(whichInput == "coast" & (k.1 == 0 | k.2 > 0))) |>
  subset(!(whichInput == "depth" & (k.1 > 0 | k.2 == 0)))

ppm.full.set <- rbind(ppm.full.set, data.frame(k.1 = 0, k.2 = 0, whichInput = "none"))
ppm.full.set$basis <- "NA"

ppm.full.set.mod <- lapply(1:nrow(ppm.full.set), \(x) {
  if (ppm.full.set$k.1[x] > 2 | ppm.full.set$k.2[x] > 2) {
    toReturn <- ppm.full.set[rep(x, 2), ]
    toReturn$basis <- c("tp", "cr")
  } else {
    toReturn <- ppm.full.set[x, ]
  }

  return (toReturn)
  }) |>
  do.call(what = "rbind")

rownames(ppm.full.set.mod) <- 1:nrow(ppm.full.set.mod)

# -------------------------------------------------------------------------

ppm.formulas <- apply(ppm.full.set.mod, 1, function(x) {
  base <- "sightings.pp ~ "

  if (x[3] == "none") {
    base <- paste(base, "1")
  } else {
    if (x[1] != "0") {
      if (x[1] == "1") {
        base <- paste0(base, "coast")
      } else if (x[1] == "2") {
        base <- paste0(base, "I(coast^2)")
      } else {
        base <- paste0(base, "s(coast, k = ", x[1], ", bs = \"", x[4], "\")")
      }
    }

    if (x[3] == "both") {
      base <- paste0(base, " + ")
    }

    if (x[2] != "0") {
      if (x[2] == "1") {
        base <- paste0(base, "depth")
      } else if (x[2] == "2") {
        base <- paste0(base, "I(depth^2)")
      } else {
        base <- paste0(base, "s(depth, k = ", x[2], ", bs = \"", x[4], "\")")
      }
    }
  }

  return (as.formula(base))
})

# -------------------------------------------------------------------------

sighting.models <- lapply(ppm.formulas, ppm, use.gam = TRUE, data = sightings.covariates)

ppm.full.set.mod$modelNum <- 1:nrow(ppm.full.set.mod)
ppm.full.set.mod$valid <- sapply(sighting.models, valid)
ppm.full.set.mod$AIC <- sapply(sighting.models, AIC)

# ppm.full.set.mod.Gres <- lapply(sighting.models, envelope, fun = Gres, correction = "bord", nsim = 200, verbose = FALSE)

results.sorted <- subset(ppm.full.set.mod, valid) |>
  {\(x){
    ix <- sort(x$AIC, index.return = TRUE)$ix
    x[ix, ]
  }}()

ppm.formulas[[results.sorted$modelNum[1]]]

# -------------------------------------------------------------------------
sightings.ppm <- ppm(sightings.pp ~ s(coast, k = 3, bs = "tp"), data = sightings.covariates, use.gam = TRUE)

## Expect that the residuals sum to zero
residuals(sightings.ppm) |>
  integral.msr()

## Visualise the statistic for assessing the suitability of the model fit under the appropriate correction
# sightings.sim <- envelope(sightings.pp, Kinhom, lambda = sightings.ppm, nsim = 100, update = TRUE, correction = "isotropic")
sightings.sim <- envelope(sightings.pp, Kinhom,
                          # lambda = predict(sightings.ppm,  type="trend"),
                          simulate = expression(rpoispp(predict(sightings.ppm, type = "trend"))),
                          nsim = 1000, correction = "translate")

## We expect, under simulation, that the variability of the `K-residual` function to capture 0 for all values of `r` to say whether we have evidence against the Poisson assumption, after accounting for spatial inhomogeneity
plot(sightings.sim)

##
Kinhom(sightings.pp, predict(sightings.ppm, locations = sightings.pp, type = "trend"), correction = "translate") |>
  plot()

## Visualise the smoothed residuals
residuals(sightings.ppm) |>
  Smooth() |>
  plot()

# -------------------------------------------------------------------------

rm(list = setdiff(ls(), c("sightings.pp", "sightings.ppm" , "sightings.covariates", "ppm.full.set.mod", "sightings.sim")))

save.image("./code/whales/ppm-Info.RData")

