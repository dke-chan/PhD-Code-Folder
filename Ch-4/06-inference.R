load("./data/whales.RData")
source("./code/whales/00-aux-funcs.R")
library(Bhat); library(sp); library(ascr); library(mgcv); library(spatstat); library(sf)
library(xtable)

# -------------------------------------------------------------------------

toTabulate <- readRDS("./code/whales/ihpp-fitted.rds")

toTabulate.models.raw <- lapply(toTabulate$modelNum, function (x) {
  readRDS(paste0("./code/whales/ihpp-model-", x, ".rds"))
})

toTabulate$deltaAIC <- toTabulate$AIC -  toTabulate$AIC[1]

##
scrDateRange <- read.table("Data/capt-toa.txt") |>
  apply(2, \(x) {as.POSIXct(ifelse(x > 0, (x - 719529) * 86400, 0), origin = "1970-01-01")}) |>
  as.numeric() |>
  {\(x){x[x > 0]}}() |>
  range() |>
  as.POSIXct(origin = "1970-01-01")

localisations.full <- read.csv("Data/obis_seamap_custom_5d33f2b4aa589_20190721_011315_points_csv.csv", stringsAsFactors = FALSE)
localisations.full$date_time_2 <- as.POSIXct(localisations.full$date_time, format = "%Y-%m-%d %H:%M:%S")
localisations.df <- subset(localisations.full, scrDateRange[1] <= date_time_2 & date_time_2 <= scrDateRange[2])

coordinates(localisations.df) <- c("longitude", "latitude")
proj4string(localisations.df) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
localisations.df <- spTransform(localisations.df, CRS("+proj=utm +zone=10 +datum=WGS84"))

localisations.num <- localisations.df@coords |>
  as.data.frame() |>
  setNames(c("x", "y")) |>
  subset(min(cali.Traps$x) <= x & x <= max(cali.Traps$x) &
           min(cali.Traps$y) <= y & y <= max(cali.Traps$y)) |>
  nrow()

control <- list(label = "logLambda", est = log(localisations.num), lower = log(1), upper = log(400))
R.A.Guazzo.Int <- exp(plkhci(control, \(x){-log(dpois(localisations.num, exp(x)))}, "logLambda")) / (localisations.num/0.43) ## Note that there were 101 calls localised across all four hydrophones

load("./code/whales/ppm-Info.RData")

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

full.set <- readRDS("./code/whales/ihpp-all.rds")

# Table 2.2/2.3 Calculations ----------------------------------------------

toTab.order <- order(toTabulate$deltaAIC)
toTabulate.rotate <- toTabulate[toTab.order, ]
toTabulate.models <- toTabulate.models.raw[toTab.order]

## Delta AIC for Table 2.2
round(toTabulate.rotate$deltaAIC, 2)

# Spline form for Table 2.2
toTabulate.rotate$basis
toTabulate.rotate$covariate.input

# Numbers for mask point size
attr(cali.Mask.Sea, "area") * 16 / 100

attr(cali.Mask.Traps, "area") * nrow(cali.Mask.Traps) / 100

# Area of the trap region in km^2
diff(range(cali.Traps$x)) * diff(range(cali.Traps$y)) / 1e+6

## D_H for Table 2.3
whatFormat <- "%.3f"
paste0(sprintf(whatFormat, toTabulate.rotate$rectEst * 100 / 48), " \\ (", sprintf(whatFormat, toTabulate.rotate$rectLwr * 100 / 48), ", ", sprintf(whatFormat, toTabulate.rotate$rectUpr * 100 / 48), ")")

# Figure 2.8 --------------------------------------------------------------
pdf("./figures/02-figure08.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", xlim = c(0, 10100), ylim = c(0, 1), xaxs = "i", axes = FALSE, xlab = "", ylab = "")
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0), at = seq(0, 10000, by = 2000),
     labels = seq(0, 10, by = 2)); mtext("Distance (km)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of Detection", side = 2, line = 1.55, font = 2)

abline(h = c(0, 1), col = "grey75")

toTabulate.models[[1]] |> show.detfn(add = TRUE, lwd = 2, col = 'black')
toTabulate.models[[2]] |> show.detfn(add = TRUE, lwd = 2, col = 'tomato')
toTabulate.models[[3]] |> show.detfn(add = TRUE, lwd = 2, col = 'steelblue')

legend("topright", lty = 1, lwd = 2, bg = "white", col = c("black", "tomato", "steelblue"), xjust = 1,
       legend = c(expression("Halfnormal "*Delta*"AIC = 0.00"), expression("Halfnormal "*Delta*"AIC = 0.72"), expression("Halfnormal "*Delta*"AIC = 1.42")))

box()
dev.off()


# Figure AA for Effects ------------------------------------------------------

?predict.ascr

newdata.df <- data.frame(
  coast = seq(min(cali.covariates$coast), 12000, length.out = 1000),
  depth = seq(-900, 0, length.out = 1000)
)

show.covariate <- function(x, covariate, use.log = TRUE, xlab = covariate, ylab = expression("log("*widehat("D")["c"]*")")) {
  which.coef <- grep(covariate, x$D.betapars, invert = TRUE)

  y.hat <- predict(x, newdata = newdata.df, use.log = use.log, set.zero = which.coef)

  if (covariate == "depth") {
    newdata.df[, covariate] <- newdata.df[, covariate] * -1
  }

  par(mai = c(0.55, 0.55, 0.1, 0.1))
  plot(y.hat ~ newdata.df[, covariate], type = "l", axes = FALSE, xlab = "", ylab = "")
  if (covariate == "depth") {
    cali.covariates[, covariate] <- cali.covariates[, covariate] * -1
  }
  rug(cali.covariates[, covariate], col = rgb(0, 0, 0, 0.25), ticksize = 0.02)
  box()
  if (covariate == "coast") {
    axis(1, tck = -0.015, at = seq(0, 12000, length.out = 7), labels = seq(0, 12, length.out = 7), cex.axis = 1.1, mgp = c(1, 0.4, 0))
  } else {
    axis(1, tck = -0.015, cex.axis = 1.1, mgp = c(1, 0.4, 0))
  }
  mtext(xlab, side = 1, line = 1.75, font = 2, cex = 0.95)
  axis(2, tck = -0.015, cex.axis = 1.1, mgp = c(1, 0.4, 0)); mtext(ylab, side = 2, line = 1.5, font = 2, cex = 0.95)
  # seq(1, length(coef.names))[-which.coef]
  # which.coef
}

pdf("./figures/02-figureAA.pdf", height = 11, width = 10)
layout(matrix(1:6, 3, 2, byrow = TRUE))

show.covariate(toTabulate.models[[1]], "depth", xlab = "Bathymetry (m)")
legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 0.00     "), cex = 1.65, bty = "n")
show.covariate(toTabulate.models[[1]], "coast", xlab = "Coastline Distance (km)")
legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 0.00    "), cex = 1.65, bty = "n")

show.covariate(toTabulate.models[[2]], "depth", xlab = "Bathymetry (m)")
legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 0.72    "), cex = 1.65, bty = "n")
show.covariate(toTabulate.models[[2]], "coast", xlab = "Coastline Distance (km)")
legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 0.72    "), cex = 1.65, bty = "n")

show.covariate(toTabulate.models[[3]], "depth", xlab = "Bathymetry (m)")
legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 1.42    "), cex = 1.65, bty = "n")
show.covariate(toTabulate.models[[3]], "coast", xlab = "Coastline Distance (km)")
legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 1.42    "), cex = 1.65, bty = "n")

dev.off()

# Figure 2.9 -------------------------------------------------------------

toPlot <- (toTabulate.rotate[1:3, c(12, 14, 15)] * 100 / 48) |>
  rbind(c(0.43, R.A.Guazzo.Int))
toPlot$model <- c("HN.1", "TH.2", "HN.3", "RA")
toPlot <- toPlot[c(4, 1, 2, 3), ]

pdf("./figures/02-figure09.pdf", height = 4, width = 9)

par(mai = c(0.55, 0.55, 0.1, 0.1))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 3), ylim = c(0.3, 0.7))

segments(0.25, R.A.Guazzo.Int, 2.75, R.A.Guazzo.Int, lty = 2, col = "grey75")

segments(seq(0.5, 2.5, length.out = 4), toPlot$rectLwr, seq(0.5, 2.5, length.out = 4), toPlot$rectUpr, lwd = 2,
         col = c("orange", "black", "tomato", "steelblue"))
segments(seq(0.45, 2.45, length.out = 4), toPlot$rectLwr, seq(0.55, 2.55, length.out = 4), toPlot$rectLwr, lwd = 2,
         col = c("orange", "black", "tomato", "steelblue"))
segments(seq(0.45, 2.45, length.out = 4), toPlot$rectUpr, seq(0.55, 2.55, length.out = 4), toPlot$rectUpr, lwd = 2,
         col = c("orange", "black", "tomato", "steelblue"))
points(x = seq(0.5, 2.5, length.out = 4), y = toPlot$rectEst, cex = 1.5, pch = 16, col = c("orange", "black", "tomato", "steelblue"))

axis(1, tck = -0.01, cex.axis = 1.1, mgp = c(1, 0.5, 0), at = seq(0.5, 2.5, length.out = 4),
     labels = c("R.A. Guazzo", "Halfnormal", "Halfnormal ", "Halfnormal"));
axis(1, tick = FALSE, cex.axis = 1.1, mgp = c(1, 1.5, 0), at = seq(0.5, 2.5, length.out = 4),
     labels = c("", expression(Delta*"AIC = 0.00"), expression(Delta*"AIC = 0.72"), expression(Delta*"AIC = 1.42")));

axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext(expression(bold("Call density p/km"^"2"*" p/h")), side = 2, line = 1.4)

box()

dev.off()

# Figure 2.10 -------------------------------------------------------------
pdf("./figures/02-figure10.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(sightings.sim, main = "", xaxs = "i", yaxs = "i", axes = FALSE, mgp = c(1.55, 0, 0),
     legendavoid = FALSE,
     legendpos = "topleft",
     xlim = range(sightings.sim$r) + c(-1, 1) * diff(range(sightings.sim$r)) * 0.02,
     ylim = range(sightings.sim$obs) + c(-0.05, 1.325) * diff(range(sightings.sim$obs)))
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0))
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0))

box()

dev.off()

# Figure 2.11 -------------------------------------------------------------

cali.Traps.X.8km <- range(cali.Traps$x) + c(-1, 1) * 9000
cali.Traps.Y.8km <- range(cali.Traps$y) + c(-1, 1) * 7750

cali.toWork <- cbind(cali.Mask.Sea, cali.covariates) |>
  subset(cali.Traps.X.8km[1] <= x & x <= cali.Traps.X.8km[2] &
           cali.Traps.Y.8km[1] <= y & y <= cali.Traps.Y.8km[2])

cali.toWork$ppm <- predict(sightings.ppm, locations = as.data.frame(cali.toWork[, 1:2]), covariates = sightings.covariates)
cali.toWork$mod1 <- predict(toTabulate.models[[1]], newdata = cali.toWork)
cali.toWork$mod2 <- predict(toTabulate.models[[2]], newdata = cali.toWork)
cali.toWork$mod3 <- predict(toTabulate.models[[3]], newdata = cali.toWork)

cali.toWork$ppm <- cali.toWork$ppm / sum(cali.toWork$ppm)
cali.toWork$mod1 <- cali.toWork$mod1 / sum(cali.toWork$mod1)
cali.toWork$mod2 <- cali.toWork$mod2 / sum(cali.toWork$mod2)
cali.toWork$mod3 <- cali.toWork$mod3 / sum(cali.toWork$mod3)

num.of.levels <- 25
cali.levels <- with(cali.toWork, cut(c(ppm, mod1, mod2, mod3), pretty(c(ppm, mod1, mod2, mod3), num.of.levels)))

cali.toWork[, c("ppm_c", "mod1_c", "mod2_c", "mod3_c")] <- as.numeric(cali.levels)

m.pixel.length <- (attr(cali.Mask.Sea, "area") * 10000) |> sqrt()
cali.toWork$xl <- cali.toWork$x - m.pixel.length/2
cali.toWork$xr <- cali.toWork$x + m.pixel.length/2
cali.toWork$yb <- cali.toWork$y - m.pixel.length/2
cali.toWork$yt <- cali.toWork$y + m.pixel.length/2

# Spectral; Vik; Lajolla; Mako; Rocket
cali.colours <- hcl.colors(num.of.levels, "Mako", rev = TRUE)
box.colour <- rgb(234, 230, 187, maxColorValue = 255)

pdf("./figures/02-figure11.pdf", height = 11.5, width = 11)
layout(matrix(c(rep(1:2, times = 4), rep(3:4, times = 4), 5, 5), nrow = 9, ncol = 2, byrow = TRUE))

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 7350,
     ylim = range(cali.Traps$y) + c(-1, 1) * 7350)
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2)

with(cali.toWork, rect(xl, yb, xr, yt, col = cali.colours[mod1_c], border = cali.colours[mod1_c]))

polygon(coastShape, lwd = 2, col = "khaki")
polygon(coastShape - matrix(rep(c(2000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5, border = "orangered")
polygon(coastShape - matrix(rep(c(4000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5)

legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 0.00"), cex = 1.75, title.adj = 0.9, bty = "n", title = "Acoustic SCR")

box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 7350,
     ylim = range(cali.Traps$y) + c(-1, 1) * 7350)
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2)

with(cali.toWork, rect(xl, yb, xr, yt, col = cali.colours[mod2_c], border = cali.colours[mod2_c]))

polygon(coastShape, lwd = 2, col = "khaki")
polygon(coastShape - matrix(rep(c(2000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5, border = "orangered")
polygon(coastShape - matrix(rep(c(4000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5)

legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 0.72"), cex = 1.75, title.adj = 0.9, bty = "n", title = "Acoustic SCR")

box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 7350,
     ylim = range(cali.Traps$y) + c(-1, 1) * 7350)
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2)

with(cali.toWork, rect(xl, yb, xr, yt, col = cali.colours[mod3_c], border = cali.colours[mod3_c]))

polygon(coastShape, lwd = 2, col = "khaki")
polygon(coastShape - matrix(rep(c(2000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5, border = "orangered")
polygon(coastShape - matrix(rep(c(4000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5)

legend("topright", legend = expression("Halfnormal "*Delta*"AIC = 1.42"), cex = 1.75, title.adj = 0.9, bty = "n", title = "Acoustic SCR")

box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 7350,
     ylim = range(cali.Traps$y) + c(-1, 1) * 7350)
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2)

with(cali.toWork, rect(xl, yb, xr, yt, col = cali.colours[ppm_c], border = cali.colours[ppm_c]))

polygon(coastShape, lwd = 2, col = "khaki")
polygon(coastShape - matrix(rep(c(2000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5, border = "orangered")
polygon(coastShape - matrix(rep(c(4000, 0), each = nrow(coastShape)), nrow = nrow(coastShape), ncol = 2), lwd = 3, col = NA, lty = 5)

legend("topright", legend = "IHPP", cex = 1.75, title.adj = 0, bty = "n", title = "Visual Survey")

box()

par(mai = c(0, 1, 0, 1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = c(0, 10),
     ylim = c(0, 1))

text(5, 0.9, "Probability Density", cex = 2, font = 2, pos = NULL)

rect.length <- (9.5 - 0.5)/num.of.levels

rect(seq(0.5, 9.5 - rect.length, length.out = num.of.levels), 0.45, seq(0.5 + rect.length, 9.5, length.out = num.of.levels), 0.725,
     col = cali.colours, border = cali.colours)
rect(0.5, 0.45, 9.5, 0.725)
segments(x0 = seq(0.5, 9.5, by = rect.length * 5), y0 = 0.35, y1 = 0.45)

# levels(cali.levels)
text(x = seq(0.5, 9.5, by = rect.length * 5), y = 0.35, pos = 1, cex = 1.75,
     labels = seq(0, 0.0125, length.out = 6))

dev.off()

# Appendix A.2 ------------------------------------------------------------

scr.full.set <- full.set
scr.full.set$AIC[which(scr.full.set$AIC == 0)] <- NA

scr.fullresults.sorted <- scr.full.set |>
  {\(x){
    ix <- sort(x$AIC, na.last = TRUE, index.return = TRUE)$ix
    x[ix, ]
  }}()

scr.nrow <- nrow(scr.fullresults.sorted)

scr.toPrint <- data.frame(DetFn = character(scr.nrow), Spline = character(scr.nrow), Bathymetry = character(scr.nrow),  Coastline = character(scr.nrow), det = numeric(scr.nrow), Std.Star = character(scr.nrow))

scr.toPrint$DetFn <- sapply(scr.fullresults.sorted$detfn, \(x) {
  # switch(x, hn = '\\texttt{"hn"}', th = '\\texttt{"th"}', hr = '\\texttt{"hr"}', hhn = '\\texttt{"hhn"}', hn.f = '\\texttt{"hn"} ($g_0 \\! = \\! 1$)')
  switch(x, hn = '\\texttt{"hn"}', th = '\\texttt{"th"}', hr = '\\texttt{"hr"}', hhn = '\\texttt{"hhn"}', hn.f = '\\texttt{"hn\\_f"}')
})

scr.temp <- round(scr.fullresults.sorted$AIC - min(scr.fullresults.sorted$AIC, na.rm = TRUE), 2)
# scr.toPrint$AIC <- ifelse(is.na(scr.temp), "---", sprintf("%.2f", scr.temp))
scr.toPrint.Included <- ifelse(scr.fullresults.sorted$noWarn == 0 & !is.nan(scr.fullresults.sorted$rectStdErr) & !is.na(scr.fullresults.sorted$AIC), "\\cmark", "\\xmark")

scr.toPrint$AIC.Star <- paste0(ifelse(is.na(scr.temp), "---", sprintf("%.2f", scr.temp)), " \\hspace{2.5mm} ", scr.toPrint.Included)

scr.toPrint$Bathymetry <- ifelse(scr.fullresults.sorted$covariate.input != "coast", paste0("$K_B = ", scr.fullresults.sorted$k.2, "$"), "---")
scr.toPrint$Coastline <- ifelse(scr.fullresults.sorted$covariate.input != "depth", paste0("$K_C = ", scr.fullresults.sorted$k.1, "$"), "---")

scr.toPrint$Spline <- ifelse(scr.fullresults.sorted$basis == "tp", "Thin-plate", "Cubic")

scr.temp <- signif(unlist(scr.fullresults.sorted$det), 2)
scr.toPrint$det <- ifelse(is.na(scr.temp), "---", sprintf(fmt = "%1.1e", scr.temp))

scr.temp <- scr.fullresults.sorted$rectStdErr |> sprintf(fmt = "%.2f")
scr.toPrint.std.err <- ifelse(scr.temp == "0.00", "---", scr.temp)

# scr.toPrint$Std.Star <- paste0(scr.toPrint.std.err, " \\hspace{2.5mm} ", scr.toPrint.Included)
scr.toPrint$Std.Star <- scr.toPrint.std.err

#
bath.linear <- which(scr.toPrint$Bathymetry == "$K_B = 1$" & scr.toPrint$Coastline == "---")
scr.toPrint$Spline[bath.linear] <- "---"
scr.toPrint <- scr.toPrint[-bath.linear[seq(2, length(bath.linear), by = 2)], ]

#
coast.linear <- which(scr.toPrint$Bathymetry == "---" & scr.toPrint$Coastline == "$K_C = 1$")
scr.toPrint$Spline[coast.linear] <- "---"
scr.toPrint <- scr.toPrint[-coast.linear[seq(2, length(coast.linear), by = 2)], ]

#
both.linear <- which(scr.toPrint$Bathymetry == "$K_B = 1$" & scr.toPrint$Coastline == "$K_C = 1$")
scr.toPrint$Spline[both.linear] <- "---"
scr.toPrint <- scr.toPrint[-both.linear[seq(2, length(both.linear), by = 2)], ]

#
bath.quad <- which(scr.toPrint$Bathymetry == "$K_B = 2$" & scr.toPrint$Coastline == "---")
scr.toPrint$Spline[bath.quad] <- "---"
scr.toPrint <- scr.toPrint[-bath.quad[seq(2, length(bath.quad), by = 2)], ]

#
coast.quad <- which(scr.toPrint$Bathymetry == "---" & scr.toPrint$Coastline == "$K_C = 2$")
scr.toPrint$Spline[coast.quad] <- "---"
scr.toPrint <- scr.toPrint[-coast.quad[seq(2, length(coast.quad), by = 2)], ]

#
both.quad <- which(scr.toPrint$Bathymetry == "$K_B = 2$" & scr.toPrint$Coastline == "$K_C = 2$")
scr.toPrint$Spline[both.quad] <- "---"
scr.toPrint <- scr.toPrint[-both.quad[seq(2, length(both.quad), by = 2)], ]

#
both.quad.two <- which(scr.toPrint$Bathymetry == "$K_B = 2$" & scr.toPrint$Coastline == "$K_C = 1$")
scr.toPrint$Spline[both.quad.two] <- "---"
scr.toPrint <- scr.toPrint[-both.quad.two[seq(2, length(both.quad.two), by = 2)], ]

#
both.quad.three <- which(scr.toPrint$Bathymetry == "$K_B = 1$" & scr.toPrint$Coastline == "$K_C = 2$")
scr.toPrint$Spline[both.quad.three] <- "---"
scr.toPrint <- scr.toPrint[-both.quad.three[seq(2, length(both.quad.three), by = 2)], ]

#
intercept <- which(scr.toPrint$Bathymetry == "$K_B = 0$" & scr.toPrint$Coastline == "$K_C = 0$")
scr.toPrint$Spline[intercept] <- "---"

scr.Table <- xtable(unique(scr.toPrint), caption = "toInsert", label = "toInsert")

print(scr.Table, sanitize.text.function = \(x){x}, booktabs = TRUE, caption.placement = "top",
      tabular.environment = "longtable", floating = FALSE, include.rownames = FALSE)

# Appendix A.3 ------------------------------------------------------------

load("./code/whales/ppm-Info.RData")

ppm.full.set.sorted <- ppm.full.set.mod |>
  {\(x){
    ix <- sort(x$AIC, index.return = TRUE)$ix
    x[ix, ]
  }}()

ppm.nrow <- nrow(ppm.full.set.sorted)

ppm.toPrint <- data.frame(Spline = character(ppm.nrow), Bathymetry = character(ppm.nrow),  Coastline = character(ppm.nrow), deltaAIC = numeric(ppm.nrow))

ppm.toPrint$deltaAIC <- round(ppm.full.set.sorted$AIC - ppm.full.set.sorted$AIC[1], 2)
ppm.toPrint$Bathymetry <- ifelse(ppm.full.set.sorted$whichInput != "coast", paste0("$K_B = ", ppm.full.set.sorted$k.2, "$"), "---")
ppm.toPrint$Coastline <- ifelse(ppm.full.set.sorted$whichInput != "depth", paste0("$K_C = ", ppm.full.set.sorted$k.1, "$"), "---")
ppm.toPrint$Spline <- ifelse(ppm.full.set.sorted$basis == "tp", "Thin-plate", ifelse(ppm.full.set.sorted$basis == "cr", "Cubic", "---"))

ppm.Table <- xtable(ppm.toPrint, caption = "toInsert", label = "toInsert", align = "llccr")

nrow(ppm.Table)

print(ppm.Table, sanitize.text.function = \(x){x}, include.rownames = FALSE, booktabs = TRUE, caption.placement = "top", tabular.environment = "longtable", floating = FALSE)

# -------------------------------------------------------------------------

# set.seed(27874848)
# n.boots <- 250
# rng.seeds <- sample(0:1e9, n.boots, FALSE)

# if (!file.exists("./code/whales/ihhp-bootstrap-flag.rds")) {
#   cl <- makeCluster(3)
#   clusterExport(cl, c("covariates.traps.df"))
#   clusterEvalQ(cl, {
#     library(ascr); library(secr)
#     source("./code/whales/00-aux-funcs.R")
#   })
#   for (i in 1:nrow(toTabulate)) {
#     i.boot.fit <- parLapplyLB(cl, rng.seeds, boostrapFit, whichSource = toTabulate.models[[i]], whichThresh = 16000)
#     saveRDS(i.boot.fit, paste0("./code/whales/ihpp-model-boots-", i, ".rds"))
#   }
#
#   stopCluster(cl)
#
#   saveRDS(0, "./code/whales/ihhp-bootstrap-flag.rds")
# } else {
#
# }
