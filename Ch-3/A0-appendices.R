library(ascr); library(dplyr)
source("./code/sonobuoys/00-aux-funcs.R")
load("./data/sonobuoys.RData")

detfn.choice <- "hhn"
obs.pointers <- which(rowSums(full.capt.hist.mat) > 0)

# -------------------------------------------------------------------------
start.mod <- fit.ascr(capt = collapsed.capt, traps = buoy.start.xy, mask = mask.start, detfn = detfn.choice)
start.mod.abs <- fit.ascr(capt = collapsed.capt, traps = buoy.abs.start.xy, mask = mask.abs.start, detfn = detfn.choice)
drift.mod <- fit.ascr(capt = capt.unbatched.mrds, traps = buoy.unbatched.xy, mask = mask.unbatched, detfn = detfn.choice)
if (detfn.choice == "hhn") {
  mask.obj.session.mod <- seaTrial.df[, c("buoy.rel.x", "buoy.rel.y")] |>
    as.matrix() |>
    create.mask(buffer = 6000, spacing = sqrt(attr(mask.obj, "area") * 10000)) |>
    list() |>
    rep(length(capt.hist.session))
  drift.mod.batch <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session, mask = mask.obj.session.mod, detfn = "hhn")
} else {
  drift.mod.batch <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session, mask = mask.obj.session, detfn = detfn.choice)
}

full.boat.ave.utm.mat <- seaTrial.df %>%
  group_by(PlaybackNumber) %>%
  summarise(boat.ave.utm.x = mean(boat.utm.x), boat.utm.rel.y = mean(boat.utm.y)) %>%
  select(-PlaybackNumber) %>%
  as.matrix()

session.pointers <- data.frame(batches = factor(group.vec$TrackedGroup), detected = rowSums(full.capt.hist.mat)) %>%
  group_by(batches) %>%
  filter(detected > 0) %>%
  mutate(NewNumber = 1:n(),
         batches = as.numeric(batches))

estimate.plottr <- function(whichOb, whichPair = "start", horizontal = TRUE) {
  if (length(whichOb) > 2) {
    stop("I can't print more than four plots!")
  }

  if (length(whichOb) > 1) {
    if (horizontal) {
      cbind(rbind(rep(1, 8), matrix(c(rep(2, 8 * 4), rep(3, 8 * 4)), 8, 8)),
            rbind(rep(4, 8), matrix(c(rep(5, 8 * 4), rep(6, 8 * 4)), 8, 8))) |>
        layout()
    } else {
      rbind(rbind(rep(1, 8), matrix(c(rep(2, 8 * 4), rep(3, 8 * 4)), 8, 8)),
            rbind(rep(4, 8), matrix(c(rep(5, 8 * 4), rep(6, 8 * 4)), 8, 8))) |>
        layout()
    }
  } else {
    layout(rbind(rep(1, 8), matrix(c(rep(2, 8 * 4), rep(3, 8 * 4)), 8, 8)))
  }

  for (i in whichOb) {
    par(mai = c(0, 0, 0, 0))
    plot(0, xlim = c(0, 10), ylim = c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
    text(5, 0.5, labels = paste0("Detected Call #", i), cex = 2, font = 2)

    par(mai = c(0.4, 0.4, 0, 0.1))
    if (whichPair == "start") {
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
           xlim = range(c(seaTrial.df$buoy.utm.x, seaTrial.df$boat.utm.x)) + 500 * c(-1, 1),
           ylim = range(c(seaTrial.df$buoy.utm.y, seaTrial.df$boat.utm.y)) + 500 * c(-1, 1))

      locations(start.mod.abs, i, add = TRUE, show.labels = FALSE, trap.col = "tomato",
                cols = list(combined = "steelblue"), lty = list(combined = "longdash"))
      points(buoy.abs.start.xy, pch = 21, bg = "tomato", cex = 1.25)

      points(full.boat.ave.utm.mat[obs.pointers[i], 1], full.boat.ave.utm.mat[obs.pointers[i], 2], pch = 21, bg = "steelblue", cex = 1.5)

      axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 11S)", side = 1, line = 1.55, font = 2, cex = 0.75)
      axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 11S)", side = 2, line = 1.55, font = 2, cex = 0.8)

      box()

      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
           xlim = range(c(seaTrial.df$buoy.rel.x, seaTrial.df$boat.rel.x)) + 500 * c(-1, 1),
           ylim = range(c(seaTrial.df$buoy.rel.y, seaTrial.df$boat.rel.y)) + 500 * c(-1, 1))

      locations(start.mod, i, add = TRUE, show.labels = FALSE, trap.col = "tomato",
                cols = list(combined = "steelblue"), lty = list(combined = "longdash"))
      points(buoy.start.xy, pch = 21, bg = "tomato", cex = 1.25)

      points(full.boat.ave.rel.mat[obs.pointers[i], 1], full.boat.ave.rel.mat[obs.pointers[i], 2], pch = 21, bg = "steelblue", cex = 1.5)

      axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (Relative)", side = 1, line = 1.55, font = 2, cex = 0.75)
      axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (Relative)", side = 2, line = 1.55, font = 2, cex = 0.8)

      box()
    } else {
      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
           xlim = range(c(seaTrial.df$buoy.rel.x, seaTrial.df$boat.rel.x)) + 500 * c(-1, 1),
           ylim = range(c(seaTrial.df$buoy.rel.y, seaTrial.df$boat.rel.y)) + 500 * c(-1, 1))

      locations(drift.mod.batch, id = session.pointers$NewNumber[i], session = session.pointers$batches[i],
                add = TRUE, show.labels = FALSE, trap.col = "tomato",
                cols = list(combined = "steelblue"), lty = list(combined = "longdash"))

      points(buoy.ave.locs.session[[session.pointers$batches[i]]], pch = 21, bg = "tomato", cex = 1.25)

      points(full.boat.ave.rel.mat[obs.pointers[i], 1], full.boat.ave.rel.mat[obs.pointers[i], 2], pch = 21, bg = "steelblue", cex = 1.5)

      axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (Relative—Batched)", side = 1, line = 1.55, font = 2, cex = 0.75)
      axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (Relative—Batched)", side = 2, line = 1.55, font = 2, cex = 0.8)

      box()

      plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
           xlim = range(c(seaTrial.df$buoy.rel.x, seaTrial.df$boat.rel.x)) + 500 * c(-1, 1),
           ylim = range(c(seaTrial.df$buoy.rel.y, seaTrial.df$boat.rel.y)) + 500 * c(-1, 1))

      locations(drift.mod, session = i, show.labels = FALSE, trap.col = "tomato",
                infotypes = "combined", cols = "steelblue", lty = "longdash", add = TRUE)

      points(buoy.unbatched.xy[[i]], pch = 21, bg = "tomato", cex = 1.25)

      points(full.boat.ave.rel.mat[obs.pointers[i], 1], full.boat.ave.rel.mat[obs.pointers[i], 2], pch = 21, bg = "steelblue", cex = 1.5)

      axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (Relative—Unbatched)", side = 1, line = 1.55, font = 2, cex = 0.75)
      axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (Relative—Unbatched)", side = 2, line = 1.55, font = 2, cex = 0.8)

      box()
    }
  }
}


# Figure 3.9  -------------------------------------------------------------
pdf("./figures/03-figure09.pdf", height = 9, width = 8)
estimate.plottr(c(1, 23), "start", FALSE)
dev.off()

# Appendix A --------------------------------------------------------------
pdf("./figures/3A/start-01.pdf", height = 4.5, width = 8)
estimate.plottr(1, "start", TRUE)
dev.off()

fileNames <- paste0("./figures/3A/start-", sprintf("%02d", 2:12), ".pdf")
for (i in 1:11) {
  pdf(fileNames[i], height = 9, width = 8)
  estimate.plottr(c(i*2, i*2 + 1), "start", FALSE)
  dev.off()
}

# Appendix B --------------------------------------------------------------
pdf("./figures/3A/drift-01.pdf", height = 4.5, width = 8)
estimate.plottr(1, "drift", TRUE)
dev.off()

fileNames <- paste0("./figures/3A/drift-", sprintf("%02d", 2:12), ".pdf")
for (i in 1:11) {
  pdf(fileNames[i], height = 9, width = 8)
  estimate.plottr(c(i*2, i*2 + 1), "drift", FALSE)
  dev.off()
}
