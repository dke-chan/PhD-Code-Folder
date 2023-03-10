library(sp); library(dplyr); library(lubridate); library(tidyr); library(ascr)

# -------------------------------------------------------------------------

raw.data <- read.csv("data/DIFAR_playback.csv", header = TRUE)

buoy.loc <- SpatialPoints(cbind(raw.data$BuoyLong, raw.data$BuoyLat), proj4string = CRS("+proj=longlat"))
buoy.utm <- spTransform(buoy.loc, CRS("+proj=utm +zone=11S"))@coords
boat.loc <- SpatialPoints(cbind(raw.data$BoatLong, raw.data$BoatLat), proj4string = CRS("+proj=longlat"))
boat.utm <- spTransform(boat.loc, CRS("+proj=utm +zone=11S"))@coords

# -------------------------------------------------------------------------

seaTrial.df <- raw.data %>%
  mutate(
    buoy.utm.x = buoy.utm[, 1],
    buoy.utm.y = buoy.utm[, 2],
    boat.utm.x = boat.utm[, 1],
    boat.utm.y = boat.utm[, 2],
    DifarFixed = if_else(!is.na(DifarFixed), DifarFixed, 0),
    DifarFixed = if_else(DifarFixed >= 360, DifarFixed - 360, DifarFixed),
    obs.bearing.rad = 2 * pi * DifarFixed / 360,
    act.bearing.rad = 2 * pi * RealBearing / 360
  ) %>%
  filter(Species == "dn1") %>%
  arrange(Buoy, PlaybackNumber) %>%
  select(-X, -UTCMilliseconds, -ClipLength, -CallType, -Intensity, -MatchedAngles, -Species, -DIFARBearing, -DifarFrequency, -SignalAmplitude, -DifarGain, -DifarFixed, -RealBearing, -BuoyLat, -BuoyLong, -BoatLat, -BoatLong, -BuoyHeading.y)

rm(raw.data, buoy.loc, buoy.utm, boat.loc, boat.utm)

# -------------------------------------------------------------------------

buoy.start.locs <- seaTrial.df %>%
  filter(PlaybackNumber == 1) %>%
  select(Buoy, buoy.utm.x, buoy.utm.y) %>%
  mutate(buoy.rel.x = buoy.utm.x - buoy.utm.x[1],
         buoy.rel.y = buoy.utm.y - buoy.utm.y[1])

NE.buoy.locs <- seaTrial.df %>%
  filter(Buoy == "NE") %>%
  select(buoy.utm.x, buoy.utm.y)

seaTrial.df <- seaTrial.df %>%
  mutate(buoy.rel.x = buoy.start.locs$buoy.rel.x[1] + buoy.utm.x - NE.buoy.locs$buoy.utm.x,
         buoy.rel.y = buoy.start.locs$buoy.rel.y[1] + buoy.utm.y - NE.buoy.locs$buoy.utm.y,
         boat.rel.x = buoy.start.locs$buoy.rel.x[1] + boat.utm.x - NE.buoy.locs$buoy.utm.x,
         boat.rel.y = buoy.start.locs$buoy.rel.y[1] + boat.utm.y - NE.buoy.locs$buoy.utm.y)

# -------------------------------------------------------------------------

full.capt.hist.mat <- seaTrial.df %>%
  select(Buoy, Detected, PlaybackNumber) %>%
  mutate(Detected = as.numeric(Detected)) %>%
  pivot_wider(names_from = Buoy, values_from = Detected) %>%
  select(-PlaybackNumber) %>%
  as.matrix()

capt.hist.mat <- full.capt.hist.mat[rowSums(full.capt.hist.mat) > 0, ]

full.bearing.mat <- seaTrial.df %>%
  select(Buoy, obs.bearing.rad, PlaybackNumber) %>%
  pivot_wider(names_from = Buoy, values_from = obs.bearing.rad) %>%
  select(-PlaybackNumber) %>%
  as.matrix()

bearing.mat <- full.bearing.mat[rowSums(full.capt.hist.mat) > 0, ]

full.distance.mat <- seaTrial.df %>%
  select(Buoy, Distance, PlaybackNumber) %>%
  pivot_wider(names_from = Buoy, values_from = Distance) %>%
  select(-PlaybackNumber) %>%
  as.matrix()

full.true.bearing.mat <- seaTrial.df %>%
  select(Buoy, act.bearing.rad, PlaybackNumber) %>%
  pivot_wider(names_from = Buoy, values_from = act.bearing.rad) %>%
  select(-PlaybackNumber) %>%
  as.matrix()

full.boat.ave.rel.mat <- seaTrial.df %>%
  group_by(PlaybackNumber) %>%
  summarise(boat.ave.rel.x = mean(boat.rel.x), boat.ave.rel.y = mean(boat.rel.y)) %>%
  select(-PlaybackNumber) %>%
  as.matrix()

boat.ave.rel.mat <- full.boat.ave.rel.mat[rowSums(full.capt.hist.mat) > 0, ]

mask.obj <- seaTrial.df %>%
  select(buoy.rel.x, buoy.rel.y) %>%
  as.matrix() %>%
  create.mask(buffer = 5000)

##
group.vec <- seaTrial.df %>%
  filter(Buoy == "NE") %>%
  select(TrackedGroup) %>%
  mutate(TrackedGroup = LETTERS[as.numeric(gsub("S([0-9]+)", "\\1", TrackedGroup))])

## (Line 73)
surveyLength.session <- seaTrial.df %>%
  select(TrackedGroup, posixDate) %>%
  mutate(TrackedGroup = LETTERS[as.numeric(gsub("S([0-9]+)", "\\1", TrackedGroup))]) %>%
  filter(complete.cases(.)) %>%
  mutate(dateTime = lubridate::parse_date_time(posixDate, "dmYHM"),
         minute = lubridate::minute(dateTime)) %>%
  group_by(TrackedGroup) %>%
  summarise(duration = max(minute) - min(minute)) %>%
  mutate(duration = ifelse(duration == 0, 1, duration)) %>%
  pull(duration)

##
capt.hist.session <- full.capt.hist.mat %>%
  as_tibble() %>%
  mutate(detected = rowSums(full.capt.hist.mat) > 0) %>%
  bind_cols(group.vec) %>%
  filter(detected) %>%
  select(-detected) %>%
  split(.$TrackedGroup) %>%
  lapply(function(session) {
    session %>%
      select(-TrackedGroup) %>%
      as.matrix()
  })

##
bearing.session <- full.bearing.mat %>%
  as_tibble() %>%
  mutate(detected = rowSums(full.capt.hist.mat) > 0) %>%
  bind_cols(group.vec) %>%
  filter(detected) %>%
  select(-detected) %>%
  split(.$TrackedGroup) %>%
  lapply(function(session) {
    session %>%
      select(-TrackedGroup) %>%
      as.matrix()
  })

##
boat.ave.locs.session <- full.boat.ave.rel.mat %>%
  as_tibble() %>%
  mutate(detected = rowSums(full.capt.hist.mat) > 0) %>%
  bind_cols(group.vec) %>%
  filter(detected) %>%
  select(-detected) %>%
  split(.$TrackedGroup) %>%
  lapply(function(session) {
    session %>%
      select(-TrackedGroup) %>%
      as.matrix()
  })

##
mask.obj.session <- lapply(vector("list", length(capt.hist.session)), function(x) mask.obj)

##
buoy.ave.locs.session <- seaTrial.df %>%
  select(Buoy, buoy.rel.x, buoy.rel.y, PlaybackNumber) %>%
  mutate(TrackedGroup = rep(pull(group.vec, TrackedGroup), 4)) %>%
  group_by(Buoy, TrackedGroup) %>%
  summarise(ave.buoy.rel.x = mean(buoy.rel.x), ave.buoy.rel.y = mean(buoy.rel.y)) %>%
  filter(TrackedGroup %in% names(capt.hist.session)) %>%
  ungroup() %>%
  split(.$TrackedGroup) %>%
  lapply(function(session) {
    session %>%
      select(ave.buoy.rel.x, ave.buoy.rel.y) %>%
      as.matrix()
  }) %>%
  unname()

# -------------------------------------------------------------------------

capt.session <- lapply(names(capt.hist.session), function (index) {
  list(
    bincapt = capt.hist.session[[index]],
    bearing = bearing.session[[index]]
  )
})

capt.session.mrds <- lapply(names(capt.hist.session), function (index) {
  list(
    bincapt = capt.hist.session[[index]],
    bearing = bearing.session[[index]],
    mrds = boat.ave.locs.session[[index]]
  )
})

capt.unbatched <- lapply(1:nrow(capt.hist.mat), \(x) {
  list(
    bincapt = capt.hist.mat[x, , drop = FALSE],
    bearing = bearing.mat[x, , drop = FALSE]
  )
})

capt.unbatched.mrds <- lapply(1:nrow(capt.hist.mat), \(x) {
  list(
    bincapt = capt.hist.mat[x, , drop = FALSE],
    bearing = bearing.mat[x, , drop = FALSE],
    mrds = full.boat.ave.rel.mat[x, , drop = FALSE]
  )
})

buoy.unbatched.xy <- lapply(which(rowSums(full.capt.hist.mat) > 0, arr.ind = TRUE), \(x) {
  seaTrial.df %>%
    select(Buoy, buoy.rel.x, buoy.rel.y, PlaybackNumber) %>%
    filter(PlaybackNumber == x) %>%
    select(buoy.rel.x, buoy.rel.y) %>%
    as.matrix()
})

mask.unbatched <- lapply(vector("list", length(capt.unbatched)), function(x) mask.obj)

# -------------------------------------------------------------------------

buoy.start.xy <- subset(seaTrial.df, PlaybackNumber == 1)[, c("buoy.rel.x", "buoy.rel.y")] |> as.matrix()
mask.start <- create.mask(buoy.start.xy, buffer = 5000)

buoy.abs.start.xy <- subset(seaTrial.df, PlaybackNumber == 1)[, c("buoy.utm.x", "buoy.utm.y")] |> as.matrix()
mask.abs.start <- create.mask(buoy.abs.start.xy, buffer = 5000)

# buoy.last.xy <- subset(seaTrial.df, PlaybackNumber == 26)[, c("buoy.rel.x", "buoy.rel.y")] |> as.matrix()
# mask.last <- create.mask(buoy.last.xy, buffer = 5000)

collapsed.capt <- list(bincapt = capt.hist.mat, bearing = bearing.mat)
collapsed.mrds <- list(bincapt = capt.hist.mat, bearing = bearing.mat, mrds = boat.ave.rel.mat)

# -------------------------------------------------------------------------

# for (i in 1:30) {
#   x <- subset(seaTrial.df, PlaybackNumber == i)
#   plot(0, type = "n", xlim = range(c(seaTrial.df$boat.utm.x, seaTrial.df$buoy.utm.x)), ylim = range(c(seaTrial.df$boat.utm.y, seaTrial.df$buoy.utm.y)))
#   points(boat.utm.y ~ boat.utm.x, data = x, col = 1:4)
#   for (j in unique(x$Buoy)) {
#     y <- subset(x, Buoy == j)
#     points(y$buoy.utm.x, y$buoy.utm.y, pch = 16)
#     # mod <- ifelse(j == 2 | j == 3, pi, 0)
#     segments(y$buoy.utm.x, y$buoy.utm.y,
#              y$buoy.utm.x + sin(y$act.bearing.rad) * y$Distance * 1000,
#              y$buoy.utm.y + cos(y$act.bearing.rad) * y$Distance * 1000)
#   }
# }
#
# rm(i, x, j, y)

# -------------------------------------------------------------------------

save.image("./data/sonobuoys.RData")
