library(ascr)

# Figure 1.1 --------------------------------------------------------------

pdf("figures/01-figure01.pdf", height = 6, width = 6)

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xaxs = "i",
     xlim = c(0, 10.1), ylim = c(0, 1))
axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of Detection", side = 2, line = 1.55, font = 2)

abline(h = c(0, 1), col = "grey75")

distance <- seq(0, 11, by = 0.1)

lines(distance, ascr:::calc.hn(distance, list(g0 = 1, sigma = 2.5)), lwd = 2, col = "tomato")
lines(distance, ascr:::calc.hn(distance, list(g0 = 0.8, sigma = 3)), lwd = 2, col = "steelblue")

legend("topright", lty = 1, lwd = 2, bg = "white", col = c("tomato", "steelblue"), xjust = 1,
       legend = c(expression("Halfnormal — g"[0]*" = 1.0, "*sigma*" = 2.5"),
                  expression("Halfnormal — g"[0]*" = 0.8, "*sigma*" = 3.0")))

box()

dev.off()

# Figure 1.2 --------------------------------------------------------------

set.seed(1234567)
N.points <- 75
N.x <- runif(N.points, -10, 10)
N.y <- runif(N.points, -10, 10)

traps.xy <- expand.grid(x = c(-1, 1), y = c(-1, 1)) |>
  as.matrix()

distance.matrix <- distances(cbind(N.x, N.y), traps.xy)
hhn.probs <- ascr:::calc.hhn(distance.matrix, list(lambda0 = 2, sigma = 2.5))
detections.matrix <- rbinom(nrow(hhn.probs) * ncol(hhn.probs), 1, as.vector(hhn.probs)) |>
  matrix(nrow = nrow(hhn.probs), ncol = ncol(hhn.probs))

detected.points <- rowSums(detections.matrix) > 0

##
survey.region <- expand.grid(x = seq(-13, 13, by = 2), y = seq(-13, 13, by = 2))
circle.radius <- 7.5
circle.pi <- seq(0, 2 * pi, length.out = 100)
c.x <- circle.radius * cos(circle.pi)
c.y <- circle.radius * sin(circle.pi)
survey.region$in.circle <- ascr::distances(as.matrix(survey.region), matrix(rep(0, 2), ncol = 2)) |>
  {\(x) x <= circle.radius}()
survey.region$xleft <- survey.region$x - 1
survey.region$xright <- survey.region$x + 1
survey.region$ybottom <- survey.region$y - 1
survey.region$ytop <- survey.region$y + 1

##
pdf("./figures/01-figure02.pdf", height = 7, width = 12)
layout(matrix(c(rep(c(1, 2), times = 9), 3, 3), nrow = 10, ncol = 2, byrow = TRUE))

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-10, 10), ylim = c(-10, 10))
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y", side = 2, line = 2, font = 2)

points(N.x, N.y, pch = 16, cex = 1.25, col = "gray75")
points(N.x[detected.points], N.y[detected.points], pch = 16, cex = 1.25, col = "steelblue")
points(traps.xy, pch = 16, cex = 2, col = "tomato")

box()

par(mai = c(0.55, 0.55, 0.1, 0.1))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(-10, 10), ylim = c(-10, 10))
axis(1, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.8, 0)); mtext("X", side = 1, line = 2.3, font = 2)
axis(2, tck = -0.01, cex.axis = 1.4, mgp = c(1, 0.7, 0)); mtext("Y", side = 2, line = 2, font = 2)

abline(v = seq(-10, 10, by = 2), h = seq(-10, 10, by = 2), lty = 5, col = "grey75")

survey.region |>
  subset(in.circle) |>
  {\(x) rect(x$xleft, x$ybottom, x$xright, x$ytop, border = "lightcoral")}()

points(N.x, N.y, pch = 16, cex = 1.25, col = "gray75")
points(N.x[detected.points], N.y[detected.points], pch = 16, cex = 1.25, col = "steelblue")
points(traps.xy, pch = 16, cex = 2, col = "tomato")

box()

par(mai = c(0.2, 1, 0, 1))
plot(0, xlim = c(0, 10), ylim = c(0, 2), type = "n", axes = FALSE, xlab = "", ylab = "")

text(2.75, 0.9, "Detector", cex = 1.5, font = 2, pos = 4)
points(x = 2.65, y = 1, cex = 2, pch = 16, col = "tomato")

text(4, 0.9, "Detected activity centre", cex = 1.5, font = 2, pos = 4)
points(x = 3.9, y = 1, cex = 2, pch = 16, col = "steelblue")

text(6.55, 0.9, "Undetected activity centre", cex = 1.5, font = 2, pos = 4)
points(x = 6.45, y = 1, cex = 2, pch = 16, col = "grey75")
# points(x = seq(1.65, 8.25, length.out = 4), y = rep(1, 4), cex = 3, pch = 15, col = rev(subset.depth.colours))
# text(x = seq(1.65, 8.25, length.out = 4), y = rep(0.9, 4), cex = 1.85, pos = 4, offset = 1,
     # labels = c("-49.6m to -65.2m", "-49.6m to -65.2m", "-49.6m to -65.2m", "-49.6m to -165.2m"))
# cut(combined.depths, num.of.cuts) |> levels()

dev.off()

