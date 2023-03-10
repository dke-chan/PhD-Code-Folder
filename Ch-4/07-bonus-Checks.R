load("./data/whales.RData")
source("./code/whales/00-aux-funcs.R")
library(Bhat); library(sp); library(ascr); library(mgcv); library(spatstat); library(sf)
library(xtable)

# -------------------------------------------------------------------------

baseShape <- st_read("./data/centralCaliShapefiles/cencal_baseline.shp") |>
  {\(x){
    sf_project("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0",
               "+proj=utm +zone=10 +datum=WGS84",
               x$geometry[[1]])
  }}()

coastShape <- rbind(
  baseShape,
  c(max(baseShape[, 1]), max(baseShape[, 2])),
  baseShape[1, ]
)

hpp.fits <- read.csv("./code/whales/hpp-fits.csv", row.names = 1) # hpp.fits <- read.csv("./code/whales/hpp-fits.csv", row.names = 1)

# -------------------------------------------------------------------------

cali.Capt <- list(bincapt = as.matrix(cali.CaptHist), toa = as.matrix(cali.ToA))
cali.covariates <- data.frame(depth = cali.Mask.z[cali.Mask.z < 0], coast = cali.Mask.minCoastDist)

tightFit.df = read.csv("./code/whales/scr-to-fit.csv")

subset(tightFit.df, detfn == "hn" & k.1 == 2 & k.2 == 4)
# k.1 = coast
# k.2 = bathy
# 67, 68, 69, 70, 72

hn.tp.23 = coreWorkerFunc(tightFit.df[216, ])
hn.tp.23$value$args$hess = TRUE
hn.tp.23 <- tryCatch(
  withWarnings(do.call("fit.ascr", hn.tp.23$value$args)),
  error = function (e) e)

cali.Traps.X.8km <- range(cali.Traps$x) + c(-1, 1) * 9000
cali.Traps.Y.8km <- range(cali.Traps$y) + c(-1, 1) * 7750

cali.toWork <- cbind(cali.Mask.Sea, cali.covariates) |>
  subset(cali.Traps.X.8km[1] <= x & x <= cali.Traps.X.8km[2] &
           cali.Traps.Y.8km[1] <= y & y <= cali.Traps.Y.8km[2])

cali.toWork$mod1 <- predict(hn.tp.23$value, newdata = cali.toWork)

num.of.levels <- 25

m.pixel.length <- (attr(cali.Mask.Sea, "area") * 10000) |> sqrt()
cali.toWork$xl <- cali.toWork$x - m.pixel.length/2
cali.toWork$xr <- cali.toWork$x + m.pixel.length/2
cali.toWork$yb <- cali.toWork$y - m.pixel.length/2
cali.toWork$yt <- cali.toWork$y + m.pixel.length/2

# Spectral; Vik; Lajolla; Mako; Rocket
cali.colours <- hcl.colors(num.of.levels, "Mako", rev = TRUE)
box.colour <- rgb(234, 230, 187, maxColorValue = 255)

# -------------------------------------------------------------------------

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 7350,
     ylim = range(cali.Traps$y) + c(-1, 1) * 7350)
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2)

heist = with(cali.toWork, cut(mod1, pretty(mod1, num.of.levels)))
with(cali.toWork, rect(xl, yb, xr, yt, col = cali.colours[heist], border = cali.colours[heist]))

polygon(coastShape, lwd = 2, col = "khaki")
polygon(coastShape - matrix(rep(c(2000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5, border = "orangered")
polygon(coastShape - matrix(rep(c(4000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5)

beepr::beep(3)
