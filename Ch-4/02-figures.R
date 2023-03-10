load("./data/whales.RData")
library(sf); library(sp); library(ascr); library(raster); library(parallel); library(spatstat); library(mgcv)

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

saveRDS(coastShape, "./code/whales/coastShape.rds")

camera.LongLat <- data.frame(latitude = 36.439892, longitude = -121.922361)
coordinates(camera.LongLat) <- c("longitude", "latitude")
proj4string(camera.LongLat) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
camera.XY <- spTransform(camera.LongLat, CRS("+proj=utm +zone=10 +datum=WGS84"))

cali.Bathymetry.LongLat <- raster("data/sccrm1s_bathymetry.tif") |>
  rasterToPoints(spatial = TRUE)
cali.Bathymetry <- (spTransform(cali.Bathymetry.LongLat, CRS("+proj=utm +zone=10 +datum=WGS84"))) |>
  {\(x){
    cbind(x@coords, z = x@data$sccrm1s_bathymetry)
  }}()
rm(cali.Bathymetry.LongLat)
gc()

# Figure 2.1 --------------------------------------------------------------

pdf("./figures/02-figure01.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 4000,
     ylim = range(cali.Traps$y) + c(-1, 1) * 4000)
# xlim = range(cali.Traps$x) + c(-1, 1) * 17500,
# ylim = range(cali.Traps$y) + c(-1, 1) * 17500)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 10)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 1.55, font = 2)

polygon(coastShape, lwd = 2, col = "khaki")
polygon(rep(range(cali.Traps$x), each = 2), c(range(cali.Traps$y), rev(range(cali.Traps$y))), lty = 2)
points(cali.Traps, pch = 21, bg = "tomato", cex = 1.5, lwd = 1)
points(camera.XY, pch = 21, bg = "#47E3FF", cex = 1.5, lwd = 1)

legend("topleft", legend = c("Hydrophone", "Research Station"),  pt.cex = 1.5, pch = 21, pt.lwd = 1, pt.bg = c("tomato", "#47E3FF"))

box()
dev.off()

# Figure 2.2 --------------------------------------------------------------

set.seed(5353)
z.x <- seq(-1.05, 1.05, by = 0.1) -> z.y
spatial.covariate <- function (x, y) {
  -0.048 * x^2 + 0.02 * x + 2 * abs(y) + -1.5 * x
}
z1 <- outer(X = z.x, Y = z.y, FUN = spatial.covariate)

sim.points <- rpoispp(function(x, y) {10 * exp(spatial.covariate(x, y))}, win = owin(c(-1.05, 1.05), c(-1.05, 1.05)))

pdf("./figures/02-figure02.pdf", height = 6.5, width = 11.5)
layout(matrix(c(rep(c(1, 2), times = 9), 3, 3), nrow = 10, ncol = 2, byrow = TRUE))
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(z.y ~ z.x, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1))
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y", side = 2, line = 2, font = 2)

image(z.x, z.y, z1, useRaster = TRUE, add = TRUE)

box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(z.y ~ z.x, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c(-1, 1), ylim = c(-1, 1))
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X", side = 1, line = 2.3, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y", side = 2, line = 2, font = 2, cex = 1.1)

points(sim.points, pch = 16, cex = 1.5)
# contour(z.x, z.y, z1, add = TRUE)

box()

par(mai = c(0, 1, 0, 1))
plot(0, xlim = c(0, 13), ylim = c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")

rect(xleft = seq(0, 11, length.out = 12) + 0.5,
     ybottom = rep(1, each = 12),
     xright = seq(1, 12, length.out = 12) + 0.5,
     ytop = rep(1.9, each = 12),
     col = hcl.colors(12, "YlOrRd", rev = TRUE),
     border = col)
rect(0.5, 1, 12.5, 1.9)

text(0, 1.35, "Z", cex = 2.5, font = 2, pos = 4, adj = 1)

# points(x = rep(seq(0.25, 9, length.out = 6), 2), y = rep(c(1.5, 0.5), each = 6), cex = 3, pch = 22, bg = hcl.colors(12, "YlOrRd", rev = TRUE))

# text(x = rep(seq(0.25, 9, length.out = 6), 2), y = rep(c(1.5, 0.5), each = 6), cex = 1.75, pos = 4, offset = 1,
#      labels = cut(z1, 12) |> levels())
toCutLevel = cut(z1, 12) |> levels()
toCutLabels = unique(c(sub("\\((.+),.*", "\\1", toCutLevel), sub("[^,]*,([^]]*)\\]", "\\1", toCutLevel)))

text(x = seq(0, 12, length.out = 13) + 0.5, y = rep(0.925, 13), cex = 1.5, pos = 1,
     labels = c(toCutLabels))
# hcl.colors(12, "YlOrRd", rev = TRUE)
#

dev.off()


# Figure 2.3 --------------------------------------------------------------

pdf("./figures/02-figure03.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 15000,
     ylim = range(cali.Traps$y) + c(-1, 1) * 15000)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 10)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 1.55, font = 2)

polygon(coastShape, lwd = 2, col = "khaki")
polygon(rep(range(cali.Traps$x), each = 2), c(range(cali.Traps$y), rev(range(cali.Traps$y))), lty = 5)
points(cali.Traps, pch = 21, bg = "tomato", cex = 1.25, lwd = 1)
points(camera.XY, pch = 21, bg = "#47E3FF", cex = 1.25, lwd = 1)
chull(cali.Mask) |>
  {\(x){lines(cali.Mask[c(x, x[1]), ], lty = 5, lwd = 2, col = "tomato")}}()

legend("topleft", legend = c("Hydrophone", "Research Station"),  pt.cex = 1.5, pch = 21, pt.lwd = 1, pt.bg = c("tomato", "#47E3FF"))

box()
dev.off()

# Figure 2.4 --------------------------------------------------------------

m.pixel.length <- (attr(cali.Mask.Sea, "area") * 10000) |> sqrt()
xRange2.4 <- range(cali.Traps$x)
yRange2.4 <- range(cali.Traps$y)

# depth <- cali.Mask.z[cali.Mask.z < 0]

cali.Mask.Subset <- cali.Mask.Sea |>
  as.data.frame() |>
  cbind(z = cali.Mask.z[cali.Mask.z < 0]) |>
  subset(xRange2.4[1] <= x & x <= xRange2.4[2] &
           yRange2.4[1] <= y & y <= yRange2.4[2])

cali.Mask.Subset$xl <- cali.Mask.Subset$x - m.pixel.length/2
cali.Mask.Subset$xr <- cali.Mask.Subset$x + m.pixel.length/2
cali.Mask.Subset$yb <- cali.Mask.Subset$y - m.pixel.length/2
cali.Mask.Subset$yt <- cali.Mask.Subset$y + m.pixel.length/2

cali.Bathmetry.Subset <- cali.Bathymetry |>
  as.data.frame() |>
  subset(min(cali.Mask.Subset$xl) <= x & x <= max(cali.Mask.Subset$xr) &
           min(cali.Mask.Subset$yb) <= y & y <= max(cali.Mask.Subset$yt))

combined.depths <- c(cali.Bathmetry.Subset$z, cali.Mask.Subset$z)

num.of.cuts <- 4
subset.depth.colours <- colorRampPalette(c("#4A1486","#084594", "#DBDDEC"))(num.of.cuts)

pdf("./figures/02-figure04.pdf", height = 7, width = 12)
layout(matrix(c(rep(c(1, 2), times = 9), 3, 3), nrow = 10, ncol = 2, byrow = TRUE))
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x),
     ylim = range(cali.Traps$y))
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2)

points(cali.Bathmetry.Subset, pch = 16, cex = 1.1,
       col = subset.depth.colours[cut(combined.depths, num.of.cuts)[1:nrow(cali.Bathmetry.Subset)] |> as.numeric()])
abline(v = unique(c(cali.Mask.Subset$xl, cali.Mask.Subset$xr)),
       h = unique(c(cali.Mask.Subset$yb, cali.Mask.Subset$yt)), lty = 5, lwd = 1.25)
box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x),
     ylim = range(cali.Traps$y))
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X (UTM Zone 10)", side = 1, line = 2.3, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 2, font = 2, cex = 1.1)

rect(cali.Mask.Subset$xl, cali.Mask.Subset$yb, cali.Mask.Subset$xr, cali.Mask.Subset$yt,
     border = subset.depth.colours[cut(combined.depths, num.of.cuts)[nrow(cali.Bathmetry.Subset) + c(1:nrow(cali.Mask.Subset))] |> as.numeric()],
     col = subset.depth.colours[cut(combined.depths, num.of.cuts)[nrow(cali.Bathmetry.Subset) + c(1:nrow(cali.Mask.Subset))] |> as.numeric()])

abline(v = unique(c(cali.Mask.Subset$xl, cali.Mask.Subset$xr)),
       h = unique(c(cali.Mask.Subset$yb, cali.Mask.Subset$yt)), lty = 5, lwd = 1.25)
box()

par(mai = c(0.2, 1, 0, 1))
plot(0, xlim = c(0, 10), ylim = c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")

text(0, 0.9, "Bathymetry", cex = 2, font = 2, pos = 4)
points(x = seq(1.65, 8.25, length.out = 4), y = rep(1, 4), cex = 3, pch = 15, col = rev(subset.depth.colours))
text(x = seq(1.65, 8.25, length.out = 4), y = rep(0.9, 4), cex = 1.85, pos = 4, offset = 1,
     labels = c("49.6m to 65.2m", "49.6m to 65.2m", "49.6m to 65.2m", "49.6m to 165.2m"))
# cut(combined.depths, num.of.cuts) |> levels()

dev.off()

# Figure 2.5 --------------------------------------------------------------

cali.Mask.Dist <- cali.Mask.Sea |>
  as.data.frame() |>
  subset(597000 <= x & x <= 598500 &
           4026200 <= y & y <= 4027500) |>
  as.matrix()

cl <- makeCluster(3)
clusterExport(cl, c("cali.Mask.Dist", "cali.Bathymetry", "cali.Mask.Land"))
clusterEvalQ(cl, {library(ascr)})
cali.Mask.Dist.Vals <- parLapplyLB(cl, 1:nrow(cali.Mask.Dist), function (rowIndex) {
  values <- distances(cali.Bathymetry[cali.Bathymetry[, 3] >= 0, 1:2], cali.Mask.Dist[rowIndex, , drop = FALSE])
  which.min(values)
}) |>
  unlist()
cali.Mask.Dist.Vals.mask <- parLapplyLB(cl, 1:nrow(cali.Mask.Dist), function (rowIndex) {
  values <- distances(cali.Mask.Land, cali.Mask.Dist[rowIndex, , drop = FALSE])
  which.min(values)
}) |>
  unlist()
stopCluster(cl)

cali.Landmass <- (cali.Bathymetry[cali.Bathymetry[, 3] >= 0, ])[cali.Mask.Dist.Vals, ]

pdf("./figures/02-figure05.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Mask.Dist[, 1]) + c(-1, 1) * 800,
     ylim = range(cali.Mask.Dist[, 2]) + c(-1, 1) * 800)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 10)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 1.55, font = 2)

polygon(coastShape, lwd = 2, col = "khaki")

segments(cali.Mask.Dist[, 1], cali.Mask.Dist[, 2], cali.Mask.Land[cali.Mask.Dist.Vals.mask, 1], cali.Mask.Land[cali.Mask.Dist.Vals.mask, 2],
         col = "deepskyblue", lwd = 2)
segments(cali.Mask.Dist[, 1], cali.Mask.Dist[, 2], cali.Landmass[, 1], cali.Landmass[, 2],
         col = "orangered", lwd = 2)

points(cali.Mask.Dist, pch = 16, cex = 1.5)
cali.Mask.Land[cali.Mask.Dist.Vals.mask, ] |>
  as.data.frame() |>
  points(cex = 1.5, pch = 21, bg = "deepskyblue")
points(cali.Landmass, pch = 21, bg = "orangered", cex = 1.5)

legend("topleft", legend = c("Mask Point (Sea)", "Bathymetry Point (Land)", "Mask Point (Land)"),
       pt.cex = 1.5, pch = c(16, 21, 21), pt.lwd = 1, pt.bg = c("black", "orangered", "deepskyblue"))

box()
dev.off()

# Figure 2.BB --------------------------------------------------------------

pdf("./figures/02-figureBB.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))

plot(I(depth*-1) ~ coast, data = cali.covariates, axes = FALSE, xlab = "", ylab = "", xlim = c(0, 20000),
     col = rgb(0,0,0,0.25), pch = 16)
axis(1, at = pretty(cali.covariates$coast), labels = c(0, 5, 10, 15, 20), tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Coastline Distance (km)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Bathymetry (m)", side = 2, line = 1.55, font = 2)


box()
dev.off()

# Figure 2.6 --------------------------------------------------------------

cali.Matrix.minCoastDist <- xtabs(cali.Mask.minCoastDist ~ ., cali.Mask.Sea)
cali.Matrix.minCoastDist[cali.Matrix.minCoastDist == 0] <- NA
cali.Matrix.depth <- xtabs(cali.Mask.z[cali.Mask.z < 0] ~ ., cali.Mask.Sea)
cali.Matrix.depth[cali.Matrix.depth == 0] <- NA

cali.depthBreaks <- seq(0, -1600, length.out = 17)
cali.minCoastDistBreaks <- seq(0, 20000, by = 1000)

pdf("./figures/02-figure06.pdf", height = 14, width = 7)

c(rep(1:2, each = 40)) |>
  matrix(20, 4, TRUE) |>
  cbind(rep(c(3, 4), each = 10)) |>
  layout()

par(mai = c(0.575, 0.575, 0.05, 0.05), cex = 1)

#
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 17500,
     ylim = range(cali.Traps$y) + c(-1, 1) * 17500)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 10)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 1.55, font = 2)

depth.colours <- colorRampPalette(c("#4A1486","#084594", "#DBDDEC"))

points(cali.Mask.Sea, pch = 15, cex = 1.25)

image(x = as.numeric(rownames(cali.Matrix.depth)),
      y = as.numeric(colnames(cali.Matrix.depth)),
      z = cali.Matrix.depth,
      col = depth.colours(length(cali.depthBreaks) - 1),
      breaks = cali.depthBreaks,
      add = TRUE, useRaster = TRUE)

polygon(coastShape, lwd = 2, col = "khaki")

box()

#
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(cali.Traps$x) + c(-1, 1) * 17500,
     ylim = range(cali.Traps$y) + c(-1, 1) * 17500)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 10)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 1.55, font = 2)

dist.colours <- colorRampPalette(c("#810f7c", "#e0ecf4"))

points(cali.Mask.Sea, pch = 15, cex = 1.25)

image(x = as.numeric(rownames(cali.Matrix.minCoastDist)),
      y = as.numeric(colnames(cali.Matrix.minCoastDist)),
      z = cali.Matrix.minCoastDist,
      col = dist.colours(length(cali.minCoastDistBreaks) - 1),
      breaks = cali.minCoastDistBreaks,
      add = TRUE, useRaster = TRUE)

polygon(coastShape, lwd = 2, col = "khaki")

box()

#
par(mai = rep(0.05, 4), cex = 1)
plot(x = rep(0.5, 22), y = -7:14, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(-7, 14))

rect(0.1, -4:11, 0.35, -3:12, col = depth.colours(length(cali.depthBreaks) - 1), border = depth.colours(length(cali.depthBreaks) - 1))
rect(0.1, -4, 0.35, 12)
text(0.5, 12.75, "Bathymetry", cex = 1.25, font = 2) # legend text

(c(0, -1*cali.depthBreaks[seq(1, 16, by = 2) + 1])) |>
  sprintf(fmt = "%.0f") |>
  paste0("m") |>
  rev() |>
  text(x = 0.41, y = seq(-4, 12, by = 2), cex = 1.1, pos = 4)
segments(x0 = 0.35, y0 = seq(-4, 12, by = 2), x1 = 0.4)

#
par(mai = rep(0.05, 4), cex = 1)
plot(x = rep(0.5, 22), y = -7:14, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(-7, 14))

rect(0.1, seq(-4, 12 - 0.8, by = 0.8), 0.35, seq(-4 + 0.8, 12, by = 0.8), col = rev(dist.colours(length(cali.minCoastDistBreaks) - 1)), border = rev(dist.colours(length(cali.minCoastDistBreaks) - 1)))
rect(0.1, -4, 0.35, 12)
text(0.5, 13.25, "Coastline\nDistance", cex = 1.25, font = 2) # legend text

(cali.minCoastDistBreaks[seq(0, 20, by = 2) + 1] / 1000) |>
  paste0("km") |>
  rev() |>
  text(x = 0.41, y =  seq(-4, 12, by = 1.6), cex = 1.1, pos = 4)
segments(x0 = 0.35, y0 = seq(-4, 12, by = 1.6), x1 = 0.4)

dev.off()

# Figure 2.7 --------------------------------------------------------------

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

coordinates(sightings.df) <- c("longitude", "latitude")
proj4string(sightings.df) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
sightings.df <- spTransform(sightings.df, CRS("+proj=utm +zone=10 +datum=WGS84"))

pdf("./figures/02-figure07.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(sightings.df@coords[, 1]) + c(-1, 1) * 1250,
     ylim = range(sightings.df@coords[, 2]) + c(-1, 1) * 1250)
# xlim = range(cali.Traps$x) + c(-1, 1) * 17500,
# ylim = range(cali.Traps$y) + c(-1, 1) * 17500)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 10)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 10)", side = 2, line = 1.55, font = 2)

polygon(coastShape, lwd = 2, col = "khaki")
points(sightings.df, pch = 16, cex = 1.5, col = rgb(0, 0, 0, 0.50))
points(cali.Traps, pch = 21, bg = "tomato", cex = 1.5, lwd = 1)
points(camera.XY, pch = 21, bg = "#47E3FF", cex = 1.5, lwd = 1)

legend("topleft", legend = c("Hydrophone", "Research Station", "Sighting Location"), pt.cex = 1.5, pch = c(21, 21, 16), pt.lwd = 1, pt.bg = c("tomato", "#47E3FF", rgb(0, 0, 0, 0.50)))

box()
dev.off()

# Appendix A.1 ------------------------------------------------------------

x <- seq(-1, 11, by = 0.1)
y <- cos(x + 1/2) + x^5 / 25000
synth.df <- data.frame(x = x, y = y, r = rnorm(length(y) * 5, sd = 0.5))
synth.df$yobs <- synth.df$y + synth.df$r

pdf("./figures/2A-figure01.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(yobs ~ x, data = synth.df,axes = FALSE, xlab = "", ylab = "", pch = 16, xlim = c(0, 10), ylim = c(-3, 6),
     col = rgb(0, 0, 0, 0.75))

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

box()
dev.off()

tp.gam <- gam(yobs ~ s(x, bs = "tp", k = 5), data = synth.df)
tp.fittedMat <- predict(tp.gam, newdata = data.frame(x = unique(x)), type = "lpmatrix")

cr.gam <- gam(yobs ~ s(x, bs = "cr", k = 5), data = synth.df)
cr.fittedMat <- predict(cr.gam, newdata = data.frame(x = unique(x)), type = "lpmatrix")

library(RColorBrewer)
lines.col <- RColorBrewer::brewer.pal(5, "Set1")

pdf("./figures/2A-figure01.pdf", height = 6, width = 12)
layout(matrix(c(1, 2), ncol = 2))
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(y = c(tp.fittedMat[, -1], cr.fittedMat[, -1]), x = rep(x, times = 8), axes = FALSE, xlab = "", ylab = "", pch = 16, xlim = c(0, 10), col = "grey75", type = "n")
lines(y = tp.fittedMat[, 2], x = unique(x), col = lines.col[1], lwd = 2)
lines(y = tp.fittedMat[, 3], x = unique(x), col = lines.col[2], lwd = 2)
lines(y = tp.fittedMat[, 4], x = unique(x), col = lines.col[3], lwd = 2)
lines(y = tp.fittedMat[, 5], x = unique(x), col = lines.col[4], lwd = 2)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(y = c(tp.fittedMat[, -1], cr.fittedMat[, -1]), x = rep(x, times = 8), axes = FALSE, xlab = "", ylab = "", pch = 16, xlim = c(0, 10), col = "grey75", type = "n")
lines(y = cr.fittedMat[, 2], x = unique(x), col = lines.col[1], lwd = 2)
lines(y = cr.fittedMat[, 3], x = unique(x), col = lines.col[2], lwd = 2)
lines(y = cr.fittedMat[, 4], x = unique(x), col = lines.col[3], lwd = 2)
lines(y = cr.fittedMat[, 5], x = unique(x), col = lines.col[4], lwd = 2)
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

box()
dev.off()
