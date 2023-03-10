load("./data/sonobuoys.RData")
source("./code/sonobuoys/00-aux-funcs.R")

library(xtable); library(tidyr); library(dplyr); library(ascr)

# Figure 3.1 --------------------------------------------------------------
set.seed(5353983)

r.x <- runif(200, -3000, 3000)
r.y <- runif(length(r.x), -2000, 2000)

d.x <- -1500
d.xt <- seq(d.x, -d.x, length.out = 25)
d.y <- 250 * c(-1, 1)

trap.locs <- cbind(d.xt, rep(d.y, each = length(d.xt)))

toy.colours <- hcl.colors(length(d.xt), "viridis")

## Stationary Stuff
dist.mat.stat <- distances(cbind(r.x, r.y), cbind(rep(d.x, 2), d.y))
capt.probs.stat <- apply(dist.mat.stat, 2, halfnormal, theta = c(1, 250), link = FALSE)
binary.capt.hist.stat <- lapply(1:length(d.xt), \(index) {
  apply(capt.probs.stat, 2, \(x) {rbinom(length(x), 1, x)})
})

## Movement Stuff
dist.mat <- distances(cbind(r.x, r.y), trap.locs)
capt.probs <- apply(dist.mat, 2, halfnormal, theta = c(1, 250), link = FALSE)
binary.capt.hist <- lapply(seq(1, length(d.xt) * 2, by = 2), \(index) {
  apply(capt.probs[, c(index, index + 1)], 2, \(x) {rbinom(length(x), 1, x)})
})

##
pdf("./figures/03-figure01.pdf", height = 4.5, width = 9)
layout(rbind(matrix(c(rep(1, 14 * 5), rep(2, 14 * 9)), 14, 14), rep(3, 14)))

par(mai = c(0.4, 0.4, 0.1, 0.1))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = -1500 + 1000 * c(-1, 1),
     ylim = 1500 * c(-1, 1))

points(r.x, r.y, col = "grey75", pch = 16, cex = 2)

for (i in 1:length(binary.capt.hist.stat)) {
  y <- binary.capt.hist.stat[[i]]
  points(r.x[rowSums(y) > 0], r.y[rowSums(y) > 0], bg = toy.colours[i], pch = 21, cex = 2)
}

points(rep(d.x, 2), d.y, pch = 16, cex = 1.5, col = "red")

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.9, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

box()

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = 2500 * c(-1, 1),
     ylim = 1500 * c(-1, 1))

points(r.x, r.y, col = "grey75", pch = 16, cex = 2)

for (i in 1:length(binary.capt.hist)) {
  y <- binary.capt.hist[[i]]
  points(r.x[rowSums(y) > 0], r.y[rowSums(y) > 0], bg = toy.colours[i], pch = 21, cex = 2)
}

arrows(d.x, d.y, d.xt[length(d.xt)], lwd = 2, length = 0.15, col = "red")
points(rep(d.x, 2), d.y, pch = 16, cex = 1.5, col = "red")

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.9, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

box()

par(mai = c(0, 0, 0, 0))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 1))

text(x = 2.05, y = 0.5, label = "Time", font = 2, pos = 4, cex = 1.5)

rect.spacing <- (2.675 + seq(0, 5, length.out = length(d.xt) + 1))
rect(rect.spacing[1:length(d.xt)], 0.5, rect.spacing[2:(length(d.xt) + 1)], 0.8, col = toy.colours)

text(x = c(mean(rect.spacing[1:2]), mean(rect.spacing[length(d.xt):(length(d.xt) + 1)])), y = 0.25, label = c("1", "25"))

dev.off()

# Figure 3.2 --------------------------------------------------------------
plot.palette <- hcl.colors(30, "viridis")

x <- seaTrial.df %>%
  group_by(PlaybackNumber) %>%
  summarise(boat.ave.x = mean(boat.utm.x), boat.ave.y = mean(boat.utm.y)) %>%
  select(boat.ave.x, boat.ave.y) |>
  as.matrix()

pdf("./figures/03-figure02.pdf", height = 7, width = 7)
layout(rbind(matrix(c(rep(1, 14)), 14, 14), rep(2, 14)))
par(mai = c(0.4, 0.4, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(c(seaTrial.df$buoy.utm.x, seaTrial.df$boat.utm.x)) + 100 * c(-1, 1),
     ylim = range(c(seaTrial.df$buoy.utm.y, seaTrial.df$boat.utm.y)) + 100 * c(-1, 1))

points(seaTrial.df$buoy.utm.x, seaTrial.df$buoy.utm.y, col = "black", pch = 4, cex = 2, lwd = 5)
points(seaTrial.df$buoy.utm.x, seaTrial.df$buoy.utm.y, col = plot.palette, pch = 4, cex = 2, lwd = 2)

points(x[, 1], x[, 2], bg = plot.palette, pch = 21, cex = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 11S)", side = 1, line = 1.9, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 11S)", side = 2, line = 1.55, font = 2)

box()

par(mai = c(0, 0, 0, 0))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 1))

text(x = 1.85, y = 0.5, label = "Time", font = 2, pos = 4, cex = 1.5)

rect.spacing <- (2.675 + seq(0, 5, length.out = 31))
rect(rect.spacing[1:30], 0.5, rect.spacing[2:31], 0.8,
     col = plot.palette)

text(x = c(mean(rect.spacing[1:2]), mean(rect.spacing[30:31])), y = 0.3, label = c("1", "30"))

dev.off()

# Figure 3.3 --------------------------------------------------------------
plot.palette <- hcl.colors(30, "viridis")

x <- seaTrial.df %>%
  group_by(PlaybackNumber) %>%
  summarise(boat.ave.x = mean(boat.rel.x), boat.ave.y = mean(boat.rel.y)) %>%
  select(boat.ave.x, boat.ave.y) |>
  as.matrix()

pdf("./figures/03-figure03.pdf", height = 7, width = 7)

layout(rbind(matrix(c(rep(1, 14)), 14, 14), rep(2, 14)))
par(mai = c(0.4, 0.4, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = range(c(seaTrial.df$buoy.rel.x, seaTrial.df$boat.rel.x)) + 100 * c(-1, 1),
     ylim = range(c(seaTrial.df$buoy.rel.y, seaTrial.df$boat.rel.y)) + 100 * c(-1, 1))

points(seaTrial.df$buoy.rel.x, seaTrial.df$buoy.rel.y, col = "black", pch = 4, cex = 2, lwd = 5)
points(seaTrial.df$buoy.rel.x, seaTrial.df$buoy.rel.y, col = plot.palette, pch = 4, cex = 2, lwd = 2)

points(x[, 1], x[, 2], bg = plot.palette, pch = 21, cex = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (Relative)", side = 1, line = 1.9, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (Relative)", side = 2, line = 1.55, font = 2)

box()

par(mai = c(0, 0, 0, 0))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 10), ylim = c(0, 1))

text(x = 1.85, y = 0.5, label = "Time", font = 2, pos = 4, cex = 1.5)

rect.spacing <- (2.675 + seq(0, 5, length.out = 31))
rect(rect.spacing[1:30], 0.5, rect.spacing[2:31], 0.8,
     col = plot.palette)

text(x = c(mean(rect.spacing[1:2]), mean(rect.spacing[30:31])), y = 0.3, label = c("1", "30"))

dev.off()

# Figure 3.4 ---------------------------------------------------------------

x <- seaTrial.df %>%
  group_by(PlaybackNumber) %>%
  summarise(boat.ave.x = mean(boat.rel.x), boat.ave.y = mean(boat.rel.y)) %>%
  select(boat.ave.x, boat.ave.y) |>
  as.matrix()

all.dist <- matrix(0, 30, 4)
buoy.char <- unique(seaTrial.df$Buoy)

for (i in 1:4) {
  y <- seaTrial.df %>%
    filter(Buoy == buoy.char[i]) %>%
    select(boat.rel.x, boat.rel.y) %>%
    as.matrix()

  for (j in 1:30) {
    all.dist[j, i] <- distances(x[j, , drop = FALSE], y[j, , drop = FALSE])
  }
}

# summary(c(all.dist))
pdf("./figures/03-figure04.pdf", height = 3.5, width = 7)
par(mai = c(0.475, 0.475, 0.1, 0.1))
hist(all.dist, breaks = 20, col = "lightblue", main = "", yaxs = "i", ylim = c(0, 42), axes = FALSE)
axis(1, tck = -0.025, mgp = c(1, 0.3, 0), cex.axis = 0.8); mtext("Euclidean distance discrepancy (m)", side = 1, line = 1.325, font = 2, cex = 1.1)
axis(2, tck = -0.025, mgp = c(1, 0.45, 0), cex.axis = 0.8); mtext("Frequency", side = 2, line = 1.35, font = 2, cex = 1.1)
box()
dev.off()

# Extra for Figure 3.4 -----------------------------------------------------
hist(all.dist, breaks = 20)
# boxplot(all.dist)
# boxplot(c(all.dist))
# boxplot.stats(all.dist)

i <- 30
seaTrial.df %>%
  filter(PlaybackNumber == i) %>%
  {plot(c(.$boat.rel.x, x[i, 1]), c(.$boat.rel.y, x[i, 2]), pch = c(rep(1, 4), 16), col = c(rep("white", 4), 4))
    text(.$boat.rel.x, .$boat.rel.y, labels = .$Buoy)}

# Figure 3.5 ---------------------------------------------------------------

x <- seaTrial.df %>%
  group_by(TrackedGroup, Buoy) %>%
  summarise(buoy.ave.rel.x = mean(buoy.rel.x), buoy.ave.rel.y = mean(buoy.rel.y)) %>%
  select(TrackedGroup, Buoy, buoy.ave.rel.x, buoy.ave.rel.y) %>%
  mutate(Tracked = factor(TrackedGroup, levels = unique(seaTrial.df$TrackedGroup)),
         Tracked = as.numeric(Tracked)) %>%
  arrange(Buoy, Tracked)

all.buoy.dist <- matrix(0, 30, 4)
buoy.char <- unique(seaTrial.df$Buoy)

for (i in 1:4) {
  y <- filter(x, Buoy == buoy.char[i])[factor(seaTrial.df$TrackedGroup[1:30], levels = unique(seaTrial.df$TrackedGroup)) |> as.numeric(), ] %>%
    ungroup() %>%
    select(buoy.ave.rel.x, buoy.ave.rel.y) %>%
    as.matrix()

  z <- seaTrial.df %>%
    filter(Buoy == buoy.char[i]) %>%
    select(buoy.rel.x, buoy.rel.y) %>%
    as.matrix()

  for (j in 1:30) {
    all.buoy.dist[j, i] <- distances(y[j, , drop = FALSE], z[j, , drop = FALSE])
  }
}

# summary(c(all.buoy.dist))
pdf("./figures/03-figure05.pdf", height = 3.5, width = 7)
par(mai = c(0.475, 0.475, 0.1, 0.1))
hist(all.buoy.dist, breaks = 20, col = "lightblue", main = "", yaxs = "i", ylim = c(0, 73), axes = FALSE)
axis(1, tck = -0.025, mgp = c(1, 0.3, 0), cex.axis = 0.8); mtext("Euclidean distance discrepancy (m)", side = 1, line = 1.325, font = 2, cex = 1.1)
axis(2, tck = -0.025, mgp = c(1, 0.45, 0), cex.axis = 0.8); mtext("Frequency", side = 2, line = 1.35, font = 2, cex = 1.1)
box()
dev.off()
