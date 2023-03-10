source("./code/misspecified/01-aux-funcs.R")

survey.designs <- lapply(1, \(x) {
  list(traps.line(2, x), traps.line(3, x), traps.line(4, x), traps.equilTriangle(x),
       traps.grid(4, x), traps.grid(9, x),
       traps.grid(25, x)
  )})[[1]]


# -------------------------------------------------------------------------

pdf("figures/05-figureD1.pdf", height = 6, width = 6)

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i",
     xlim = c(0, 10.1), ylim = c(0, 1))
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of Detection", side = 2, line = 1.55, font = 2)

abline(h = c(0, 1), col = "grey75")

distance <- seq(0, 11, by = 0.1)

lines(distance, ascr:::calc.hn(distance, list(g0 = 1, sigma = 2.5)), lwd = 2, col = "black")

box()

dev.off()


# Figure 5.1 --------------------------------------------------------------

pdf("./figures/05-figure01.pdf", height = 9, width = 9)
par(cex = 1)
layout(matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE))

plot.bounds <- data.frame(
  xlwr = c(-2.5, -2, -1.5, -2.5, -2.5, -2, -1),
  xupr = c(3.5, 4, 4.5, 3.5, 3.5, 4, 5),
  ylwr = c(-3, -3, -3, -survey.designs[[4]][2, "y"] / 2 - 3, -2.5, -2, -1),
  yupr = c(3, 3, 3, survey.designs[[4]][2, "y"] / 2 + 3, 3.5, 4, 5)
) |> as.matrix()

for (i in 1:length(survey.designs)) {
  par(mai = c(0.05, 0.05, 0.05, 0.05), xaxs = "i")
  plot(survey.designs[[i]], xlim = plot.bounds[i, 1:2], ylim = plot.bounds[i, 3:4], axes = FALSE, xlab = "", ylab = "", type = "n")

  abline(h = -5:5, v = -5:5,
         col = "grey70", lty = 2)

  points(survey.designs[[i]], pch = 16, cex = 1.5)

  box()
}

dev.off()


# New figure... -----------------------------------------------------------

library(ascr)
x.dist <- 0:3000

pdf("./figures/05-figureG1.pdf", height = 12, width = 12)
palette("Dark 2")
layout(matrix(1:4, 2, 2, byrow = TRUE))

par(mai = c(0.55, 0.55, 0.15, 0.15), xaxs = "i", oma = rep(0.3, 4))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "identity", orientation = 0), col = "black", lty = 5)
lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 1, sigma = 900)), col = 4, lwd = 2)

legend("topright", legend =  "Halfnormal", bty = "n", cex = 2, text.font = 2)

box()


plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "identity", orientation = 0), col = "black", lty = 5)
lines(x.dist, ascr:::calc.hr(x.dist, list(g0 = 1, sigma = 1400, z = 8)), col = 6, lwd = 2)

legend("topright", legend =  "Hazard Rate", bty = "n", cex = 2, text.font = 2)

box()

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "identity", orientation = 0), col = "black", lty = 5)
lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 25, sigma = 550)), col = 3, lwd = 2)

legend("topright", legend =  "Hazard Halfnormal", bty = "n", cex = 2, text.font = 2)

box()

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "identity", orientation = 0), lwd = 2, col = 1)
# lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "identity", orientation = 0), col = "black", lty = 5)

box()

legend("topright", legend = "Signal Strength", bty = "n", cex = 2, text.font = 2)

dev.off()


# Figure 5.2 --------------------------------------------------------------

library(ascr)
x.dist <- 0:3000

pdf("./figures/05-figure02.pdf", height = 6, width = 12)
layout(matrix(c(1, 2), 1, 2))

par(mai = c(0.55, 0.55, 0.15, 0.15), xaxs = "i", oma = rep(0.3, 4))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 0.5, sigma = 500)), col = "tomato")
lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 0.75, sigma = 500)), col = "steelblue")
lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 1, sigma = 500)), col = "forestgreen")

legend("topright", expression(italic(g[0])), bty = "n", cex = 2)

box()

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 1, sigma = 250)), col = "tomato")
lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 1, sigma = 500)), col = "steelblue")
lines(x.dist, ascr:::calc.hn(x.dist, list(g0 = 1, sigma = 750)), col = "forestgreen")

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

legend("topright", expression(sigma), bty = "n", cex = 2)

box()

dev.off()

# Figure 5.3 --------------------------------------------------------------

pdf("./figures/05-figure03.pdf", height = 6, width = 12)
layout(matrix(c(1, 2), 1, 2))

par(mai = c(0.55, 0.55, 0.15, 0.15), xaxs = "i", oma = rep(0.3, 4))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 1.5, sigma = 500)), col = "tomato")
lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 3.5, sigma = 500)), col = "steelblue")
lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 5.5, sigma = 500)), col = "forestgreen")

legend("topright", expression(lambda[0]), bty = "n", cex = 2)

box()


plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 5.5, sigma = 250)), col = "tomato")
lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 5.5, sigma = 500)), col = "steelblue")
lines(x.dist, ascr:::calc.hhn(x.dist, list(lambda0 = 5.5, sigma = 750)), col = "forestgreen")

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

legend("topright", expression(sigma), bty = "n", cex = 2)

box()

dev.off()

# -------------------------------------------------------------------------

pdf("./figures/05-figureAA.pdf", height = 12, width = 12)
layout(matrix(1:4, 2, 2, byrow = TRUE))

par(mai = c(0.55, 0.55, 0.15, 0.15), xaxs = "i", oma = rep(0.3, 4))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1000, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue", lty = 5)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1000, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato", lty = 5)

legend("topright", expression(beta[0]), bty = "n", cex = 2)

box()


plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 1., sigma.ss = 100, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue", lty = 5)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 1., sigma.ss = 100, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato", lty = 5)

legend("topright", expression(beta[1]), bty = "n", cex = 2)

box()

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 200, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue", lty = 5)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 200, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato", lty = 5)

legend("topright", expression(sigma[s]), bty = "n", cex = 2)

box()

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(x.dist), ylim = c(0, 1))

abline(h = seq(0, 1, 0.2), v = seq(0, 3000, 500), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "spherical", orientation = 0), col = "steelblue")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "spherical", orientation = 0), col = "steelblue", lty = 5)

lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 750), ss.link = "identity", orientation = 0), col = "tomato")
lines(x.dist, ascr:::calc.ss(x.dist, list(b0.ss = 1250, b1.ss = 0.5, sigma.ss = 100, cutoff = 500), ss.link = "identity", orientation = 0), col = "tomato", lty = 5)

box()

legend("topright", expression(italic(h)), bty = "n", cex = 2)

dev.off()
