# Figure 3.1* -------------------------------------------------------------

detectors <- data.frame(x = seq(-2, 2, by = 2),  y = 0)

animal.t <- data.frame(x = c(-3, 1.3), y = c(2, -2))

animal.s <- data.frame(
  x = c(-3.5, -2.3, -4, 1.2, 1.75, 0, 3),
  y = c(0.8, 2.4, 2.75, -3, -0.8, -1.8, -2.1),
  id = c(rep(1, 3), rep(2, 4)),
  detect = c(1, 1, 0, 0, 1, 1, 0),
  arr.1 = c(1, 1, 0, 0, 0, 0, 0),
  arr.2 = c(0, 1, 0, 0, 1, 1, 0),
  arr.3 = c(0, 0, 0, 0, 1, 0, 0)
)

pdf("./figures/04-figureAA.pdf", height = 6, width = 6)
par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
     xlim = c(-4, 4), ylim = c(-4, 4))
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y", side = 2, line = 1.55, font = 2)

segments(x0 = subset(animal.s, arr.1 > 0)$x, y0 = subset(animal.s, arr.1 > 0)$y,
         x1 = detectors$x[1], y1 = detectors$y[1], lty = 5, col = "black")

segments(x0 = subset(animal.s, arr.2 > 0)$x, y0 = subset(animal.s, arr.2 > 0)$y,
         x1 = detectors$x[2], y1 = detectors$y[2], lty = 5, col = "black")

segments(x0 = subset(animal.s, arr.3 > 0)$x, y0 = subset(animal.s, arr.3 > 0)$y,
         x1 = detectors$x[3], y1 = detectors$y[3], lty = 5, col = "black")

points(detectors, cex = 1.25, pch = 16, col = "tomato")
points(animal.t, cex = 1.25, pch = 1, col = "black")
points(animal.s, cex = 1.25, pch = 16, col = ifelse(animal.s$detect == 1, "steelblue", "grey50"))

text(-2.8, 1.8, expression(bold(t)[1]))
text(-2.8, 1.9, expression(phantom(bold(t))*"*"))
text(1.5, -2.2, expression(bold(t)[2]))
text(1.5, -2.1, expression(phantom(bold(t))*"*"))

text(-1.8, -0.2, expression(bold(x)[1]))
text(-1.8, -0.1, expression(phantom(bold(x))*"*"))
text(0.2, 0.2, expression(bold(x)[2]))
text(0.2, 0.3, expression(phantom(bold(x))*"*"))
text(2.2, -0.2, expression(bold(x)[3]))
text(2.2, -0.1, expression(phantom(bold(x))*"*"))

text(-3.2, 1, expression(bold(s)[11]))
text(-3.25, 1.075, expression(phantom(bold(s))*"*"))
text(-2, 2.6, expression(bold(s)[12]))
text(-2.05, 2.675, expression(phantom(bold(s))*"*"))
text(-3.7, 2.6, expression(bold(s)[13]))
text(-3.75, 2.675, expression(phantom(bold(s))*"*"))

text(1.4, -3.3, expression(bold(s)[21]))
text(1.35, -3.225, expression(phantom(bold(s))*"*"))
text(1.95, -1, expression(bold(s)[22]))
text(1.9, -0.925, expression(phantom(bold(s))*"*"))
text(0.2, -2, expression(bold(s)[23]))
text(0.15, -1.925, expression(phantom(bold(s))*"*"))
text(3.2, -2.3, expression(bold(s)[24]))
text(3.15, -2.225, expression(phantom(bold(s))*"*"))

box()
dev.off()
