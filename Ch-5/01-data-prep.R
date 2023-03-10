library(readxl); library(ascr); library(secr); library(Rcpp); library(RcppArmadillo)

# -------------------------------------------------------------------------

load("./data/KP_ascr_dat(2021-11-29)_unfiltereddetections(forBenStevenson).RData")

# -------------------------------------------------------------------------

posts.df <- lapply(posts, \(x) {
  y <- as.data.frame(x)
  y$post <- seq(1, nrow(y))

  y
}) |>
  do.call(what = "rbind")
posts.df$array <- gsub("(0-9)*\\.[0-9]*", "\\1", rownames(posts.df)) |> as.integer()
rownames(posts.df) <- seq(1, nrow(posts.df))

detections.df <- dplyr::left_join(detections, posts.df, by = c("array", "post")) |>
  as.data.frame() |>
  unique()

# -------------------------------------------------------------------------

posts.full <- lapply(unique(posts.df$array), \(x) {
  y <- posts.df[which(posts.df$array == x), ]
  rownames(y) <- 1:nrow(y)
  return (as.matrix(y[, c("x", "y")]))
})

capt.full <- vector("list", length(posts.full))

capt.occ.1 <- create.capt(captures = subset(detections.df[, c("array", "group", "occasion", "post", "bearing", "dist")], occasion == 1),
                          traps = posts.full)
# Ignoring that individual 25 was detected by some traps more than once.
# Ignoring that individual 29 was detected by some traps more than once.

capt.occ.2 <- create.capt(captures = subset(detections.df[, c("array", "group", "occasion", "post", "bearing", "dist")], occasion == 2),
                          traps = posts.full)
# Ignoring that individual 31 was detected by some traps more than once.
# Ignoring that individual 51 was detected by some traps more than once.
# Ignoring that individual 53 was detected by some traps more than once.
# Ignoring that individual 50 was detected by some traps more than once.
# Ignoring that individual 45 was detected by some traps more than once.

capt.occ.3 <- create.capt(captures = subset(detections.df[, c("array", "group", "occasion", "post", "bearing", "dist")], occasion == 3),
                          traps = posts.full)
# Ignoring that individual 38 was detected by some traps more than once.
# Ignoring that individual 49 was detected by some traps more than once.
# Ignoring that individual 45 was detected by some traps more than once.

# Investigating multiple-detections ---------------------------------------
raw.visualiser <- function(occasion, group, whichObs = NULL) {
  x <- detections.df[which(detections.df$occasion == occasion & detections.df$group == group), ]
  theArray <- posts.df[which(posts.df$array == unique(x$array)), ]

  if (is.null(whichObs)) {
    whichObs <- seq(1, nrow(x))
  }

  plot(0, mgp = c(1.75, 0.5, 0), xlim = range(theArray$x) + 1500 * c(-1, 1), ylim = range(theArray$y) + 1500 * c(-1, 1), type = "n", xlab = "X", ylab = "Y", main = paste("Group No.", group, "Occasion No.", occasion))

  cat("Array ", unique(x$array), "\n", sep = "")
  for (i in whichObs) {
    if (i >= 1 & i <= nrow(x)) {
      arrows(x$x[i], x$y[i],
             x$x[i] + x$dist[i] * sin(x$bearing[i]),
             x$y[i] + x$dist[i] * cos(x$bearing[i]),
             lwd = 2, col = 1 + i)
      cat("Bearing ", i, " is: ", x$bearing[i], "\n", sep = "")
    }
  }

  points(theArray, pch = 16)
  legend("topright", legend = paste("Detect No.", whichObs), col = 1 + whichObs, lty = 1)
}

# 04-figure-A1 ------------------------------------------------------------

pdf("figures/04-figureA1.pdf", width = 10, height = 15)
layout(matrix(1:6, 3, 2))
par(mai = c(0.6, 0.6, 0.35, 0.15), cex = 1)
raw.visualiser(1, 25, c(1, 2, 3)) ## Keep the first two bearings (out of three bearings)


raw.visualiser(2, 31, c(1, 2, 3)) ## Keep the first and third bearings (out of three bearings)
raw.visualiser(2, 50, c(1, 2, 3)) ## Keep the first two bearings (out of three bearings)

raw.visualiser(1, 29, c(1, 2, 3, 4)) ## Keep the second and third bearings (out of four bearings)
capt.occ.1[[4]]$bearing["29", 2] <- subset(detections.df, occasion == 1 & group == 29)$bearing[2]
raw.visualiser(2, 45, c(1, 2)) ## Keep the first bearing (out of two bearings)


raw.visualiser(2, 51, c(1, 2, 3)) ## Keep the first two bearings (out of three bearings)
dev.off()

# 04-figure-A2 ------------------------------------------------------------

pdf("figures/04-figureA2.pdf", width = 10, height = 10)
layout(matrix(1:4, 2, 2))
par(mai = c(0.6, 0.6, 0.35, 0.15), cex = 1)

raw.visualiser(3, 38, c(1, 2)) ## Keep the first bearing (out of two bearings)
raw.visualiser(3, 49, c(1, 2, 3)) ## Keep the first and third bearing (out of three bearings)
raw.visualiser(3, 45, c(1, 2)) ## Keep the first bearing (out of two bearings)

raw.visualiser(2, 53, c(1, 2, 3)) ## Keep the second and third  bearings (out of three bearings)
capt.occ.2[[7]]$bearing["53", 3] <- subset(detections.df, occasion == 2 & group == 53)$bearing[2]


dev.off()

# -------------------------------------------------------------------------

capt.full <- vector("list", length(posts.full))

for (x in 1:length(posts.full)) {
  capt.full[[x]]$bincapt <- rbind(
    capt.occ.1[[x]]$bincapt, capt.occ.2[[x]]$bincapt, capt.occ.3[[x]]$bincapt
  )

  capt.full[[x]]$bearing <- rbind(
    capt.occ.1[[x]]$bearing, capt.occ.2[[x]]$bearing, capt.occ.3[[x]]$bearing
  )

  capt.full[[x]]$dist <- rbind(
    capt.occ.1[[x]]$dist, capt.occ.2[[x]]$dist, capt.occ.3[[x]]$dist
  )
}

capt.bearings <- lapply(capt.full, \(x) {x$dist <- NULL; return(x)})

raw.animal.ids <- lapply(capt.bearings, \(x) {
  rownames(x$bincapt) |> unique()
})

# secr to guess initial buffer size ---------------------------------------

secr.hists <- lapply(unique(posts.df$array), \(x) {
  y <- try(convert.capt.to.secr(capt.bearings[[x]], posts.full[[x]]))
  if (!inherits(y, "try-error")) {
    attr(y, "session") <- as.character(x)
  }

  y
})

secr.buffers <- lapply(secr.hists, \(x) {
  y <- try(suggest.buffer(x), silent = TRUE)

  if (!inherits(y, "try-error")) {
    return (y)
  }

  return (NA)
})

buffer.choice <- max(unlist(secr.buffers), na.rm = TRUE)
buffer.choice

# fit.ascr() to benchmark range of buffers --------------------------------

buffers <- seq(1, 6, by = 0.5) * buffer.choice

if (!file.exists("./code/gibbons/01-data-prep-buffer-explore.rds")) {
  buffer.fits <- lapply(buffers, \(x) {
    buffer.mask <- create.mask(posts.full, buffer = x)
    buffer.fit <- fit.ascr(capt.bearings, traps = posts.full, detfn = "hhn", mask = buffer.mask)
    return (buffer.fit)
  })

  estimated.betas <- lapply(buffer.fits, coef, pars = c("D", "lambda0", "sigma", "kappa", "esa")) |>
    do.call(what = "rbind") |>
    cbind(buffers)

  saveRDS(estimated.betas, "./code/gibbons/01-data-prep-buffer-explore.rds")
} else {
  estimated.betas <- readRDS("./code/gibbons/01-data-prep-buffer-explore.rds")
}

layout(matrix(1:4, 2, 2))
plot(D ~ buffers, data = estimated.betas, type = "l")
plot(lambda0 ~ buffers, data = estimated.betas, type = "l")
plot(sigma ~ buffers, data = estimated.betas, type = "l")
plot(kappa ~ buffers, data = estimated.betas, type = "l")
layout(1)

# Prep files to load in for analysis etc. ---------------------------------

survey.mask <- create.mask(posts.full, buffer = 5000)
bincapt.session <- lapply(capt.bearings, \(x) x$bincapt)
survey.session <- create.mask(posts.full, buffer = 5000) #, spacing = 100)
mask.dist.session <- lapply(1:length(survey.session), \(x) {ascr::distances(survey.session[[x]], posts.full[[x]])})
bearing.session <- lapply(capt.bearings, \(x) x$bearing)
expected.bearing.session <- lapply(1:length(survey.session), \(x) {t(ascr::bearings(posts.full[[x]], survey.session[[x]]))})
pixel.size.session <- sapply(survey.session, attr, "area")
animal.id.session <- lapply(capt.bearings, \(x) {
  y <- rownames(x$bincapt) |> factor() |> as.numeric()
  y - 1
})
survey.length <- 3

# Figure 4.1 --------------------------------------------------------------

pdf("./figures/04-figure01.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(posts.df$x) + c(-1, 1),
     ylim = range(posts.df$y) + c(-1, 1))
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

lapply(survey.mask, \(x) {
  chull(x) |>
    {\(y){polygon(x[c(y, y[1]), ], col = rgb(0, 0, 0, 0.1), border = NA)}}()
})

points(posts.df, pch = 21, bg = "tomato", cex = 0.9, lwd = 1)

box()
dev.off()

# -------------------------------------------------------------------------

rm(list = ls()[
  -which(ls() %in% c("survey.mask", "posts.df", "capt.bearings", "posts.full", "survey.session", "bincapt.session",
                     "survey.session", "mask.dist.session", "bearing.session", "expected.bearing.session",
                     "pixel.size.session", "animal.id.session", "survey.length", "capt.full"))
])

save.image("code/gibbons/01-data-prep.RData")
