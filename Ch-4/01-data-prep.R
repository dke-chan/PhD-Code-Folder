library(rgdal); library(sp); library(readxl); library(ascr); library(raster); library(parallel)

cali.CaptHist <- read.table("data/capt-bin.txt", col.names = c("NE", "NW", "SE", "SW"))

# The following changes the parameterisation of the ToA data to be in terms of relative ToA, as opposed to ToA since 01 November 2014
relative.ToA <- FALSE
cali.ToA <- 86400 * apply(read.table("data/capt-toa.txt"), 2, function (x) {
  x <- ifelse(x > 0, as.POSIXct(x * 86400, origin = "0000-01-01") - as.POSIXct("2014-11-01", origin = "0000-01-01"), 0)

  if (relative.ToA) {
    x.sign <- as.numeric(sign(x))
    x <- x - rep(min(x[as.logical(x.sign)])) * x.sign
  }

  return (x)}) |>
  setNames(c("NE", "NW", "SE", "SW"))

cali.Traps.LongLat <- read.table("data/trap-locs.txt", col.names = c("x", "y"))
coordinates(cali.Traps.LongLat) <- c("x", "y")
proj4string(cali.Traps.LongLat) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
cali.Traps <- as.data.frame(spTransform(cali.Traps.LongLat, CRS("+proj=utm +zone=10 +datum=WGS84")))

# cali.Bathymetry.LongLat <- read_xlsx("data/bathymetry_GraniteCanyon.xlsx", col_names = c("x", "y", "z"))
cali.Bathymetry.LongLat <- raster("data/sccrm1s_bathymetry.tif") |>
  rasterToPoints(spatial = TRUE)
# coordinates(cali.Bathymetry.LongLat) <- c("x", "y")
# proj4string(cali.Bathymetry.LongLat) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
cali.Bathymetry <- (spTransform(cali.Bathymetry.LongLat, CRS("+proj=utm +zone=10 +datum=WGS84"))) |>
  {\(x){
    cbind(x@coords, z = x@data$sccrm1s_bathymetry)
  }}()

cali.Mask.Base <- create.mask(cali.Traps, 16000)
m.pixel.length <- (attr(cali.Mask.Base, "area") * 10000) |> sqrt()
x.bound <- (min(cali.Bathymetry[, 1]) + m.pixel.length/2) <= cali.Mask.Base[, 1] &
  cali.Mask.Base[, 1] <= (max(cali.Bathymetry[, 1]) - m.pixel.length/2)
y.bound <- (min(cali.Bathymetry[, 2]) + m.pixel.length/2) <= cali.Mask.Base[, 2] &
  cali.Mask.Base[, 2] <= (max(cali.Bathymetry[, 2]) - m.pixel.length/2)
cali.Mask <- cali.Mask.Base[x.bound & y.bound, ]
attr(cali.Mask, "area") <- attr(cali.Mask.Base, "area")
attr(cali.Mask, "buffer") <- attr(cali.Mask.Base, "buffer")

cl <- makeCluster(3)
clusterExport(cl, c("cali.Mask", "m.pixel.length", "cali.Bathymetry"))
cali.Mask.z <- parLapplyLB(cl, 1:nrow(cali.Mask), function (rowIndex) {
  x.bound <- cali.Bathymetry[, 1] <= (cali.Mask[rowIndex, 1] + m.pixel.length/2) &
    (cali.Mask[rowIndex, 1] - m.pixel.length/2) <= cali.Bathymetry[, 1]
  y.bound <- cali.Bathymetry[, 2] <= (cali.Mask[rowIndex, 2] + m.pixel.length/2) &
    (cali.Mask[rowIndex, 2] - m.pixel.length/2) <= cali.Bathymetry[, 2]

  mean(cali.Bathymetry[x.bound & y.bound, 3])
}) |> unlist()
stopCluster(cl)

cali.Mask.Sea <- cali.Mask[cali.Mask.z < 0, ]
attr(cali.Mask.Sea, "area") <- attr(cali.Mask, "area")
attr(cali.Mask.Sea, "buffer") <- attr(cali.Mask, "buffer")

cali.Mask.Land <- cali.Mask[cali.Mask.z >= 0, ]
attr(cali.Mask.Land, "area") <- attr(cali.Mask, "area")
attr(cali.Mask.Land, "buffer") <- attr(cali.Mask, "buffer")

cl <- makeCluster(3)
clusterExport(cl, c("cali.Mask.Sea", "cali.Bathymetry"))
clusterEvalQ(cl, {library(ascr)})
cali.Mask.minCoastDist <- parLapplyLB(cl, 1:nrow(cali.Mask.Sea), function (rowIndex) {
  min(distances(cali.Bathymetry[cali.Bathymetry[, 3] >= 0, 1:2], cali.Mask.Sea[rowIndex, , drop = FALSE]))}) |>
  unlist()
stopCluster(cl)

cali.covariates <- data.frame(depth = cali.Mask.z[cali.Mask.z < 0], coast = cali.Mask.minCoastDist)

# Setup the components for estimating the density within the hydrophone
pixelEdgeNum <- 10
rect.xLength <- diff(range(cali.Traps$x)) / (pixelEdgeNum)
rect.yLength <- diff(range(cali.Traps$y)) / (pixelEdgeNum)

cali.Mask.Traps <- expand.grid(
  x = seq(min(cali.Traps$x) + rect.xLength/2, max(cali.Traps$x) - rect.xLength/2, length.out = pixelEdgeNum),
  y = seq(min(cali.Traps$y) + rect.yLength/2, max(cali.Traps$y) - rect.yLength/2, length.out = pixelEdgeNum)) |>
  as.matrix()
attr(cali.Mask.Traps, "area") <- rect.xLength * rect.yLength / 10000

cl <- makeCluster(3)
clusterExport(cl, c("cali.Mask.Traps", "rect.xLength", "rect.yLength", "cali.Bathymetry"))
cali.Mask.Traps.z <- parLapplyLB(cl, 1:nrow(cali.Mask.Traps), function (rowIndex) {
  x.bound <- cali.Bathymetry[, 1] <= (cali.Mask.Traps[rowIndex, 1] + rect.xLength/2) &
    (cali.Mask.Traps[rowIndex, 1] - rect.xLength/2) <= cali.Bathymetry[, 1]
  y.bound <- cali.Bathymetry[, 2] <= (cali.Mask.Traps[rowIndex, 2] + rect.yLength/2) &
    (cali.Mask.Traps[rowIndex, 2] - rect.yLength/2) <= cali.Bathymetry[, 2]

  mean(cali.Bathymetry[x.bound & y.bound, 3])
}) |> unlist()
stopCluster(cl)

cl <- makeCluster(3)
clusterExport(cl, c("cali.Mask.Traps", "cali.Bathymetry"))
clusterEvalQ(cl, {library(ascr)})
cali.Mask.Traps.minCoastDist <- parLapplyLB(cl, 1:nrow(cali.Mask.Traps), function (rowIndex) {
  min(distances(cali.Bathymetry[cali.Bathymetry[, 3] >= 0, 1:2], cali.Mask.Traps[rowIndex, , drop = FALSE]))}) |>
  unlist()
stopCluster(cl)

covariates.traps.df <- cbind(cali.Mask.Traps, depth = cali.Mask.Traps.z, coast = cali.Mask.Traps.minCoastDist) |>
  as.data.frame()

# Clean up prior to saving the data setup
remove(relative.ToA)
remove(cali.Traps.LongLat)
remove(cali.Bathymetry.LongLat)
remove(cali.Bathymetry)
remove(cali.Mask.Base)
remove(cl)
remove(x.bound)
remove(y.bound)
remove(m.pixel.length)

# Save!
save.image("./data/whales.RData")
