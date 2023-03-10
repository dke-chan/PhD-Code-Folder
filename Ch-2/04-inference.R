library(ascr); library(secr); library(xtable)

# -------------------------------------------------------------------------

files.grid <- expand.grid(scenario = "fascrNew_", design = 1:7, stringsAsFactors = FALSE)
files <- paste0(files.grid$scenario, files.grid$design, ".csv")

masterList <- data.frame(theta = numeric(0), theta.name = character(0), survey = character(0) , detfn = character(0), seed = numeric(0), aux = logical(0),
                         design = character(0))

for (index in files) {
  simData_i <- read.csv(paste0("./code/miss-v2/results/", index))
  simData_i$design <- index

  masterList <- rbind(masterList, simData_i)
}

# -------------------------------------------------------------------------

source("./code/miss-v2/01-aux-funcs.R")
sim.theta = readRDS(file = "./code/miss-v2/model_fit.RDS")

# -------------------------------------------------------------------------

halfnormal <- function(d, theta, link = TRUE) {
  if (link) {
    g0 <- exp(theta[1])/(exp(theta[1]) + 1); sigma <- exp(theta[2])
  } else {
    g0 <- theta[1]; sigma <- theta[2]
  }

  return (g0 * exp(- d^2 / (2 * sigma^2)))
}

hazardHalfnormal <- function (d, theta, link = TRUE) {
  if (link) {
    lambda0 <- exp(theta[1]); sigma <- exp(theta[2])
  } else {
    lambda0 <- theta[1]; sigma <- theta[2]
  }
  return (1 - exp(-lambda0 * exp(-(d^2 / (2*sigma^2)))))
}

hazardRate <- function (d, theta, link = TRUE) {
  if (link) {
    g0 <- exp(theta[1])/(exp(theta[1]) + 1); sigma <- exp(theta[2]); z <- exp(theta[3])
  } else {
    g0 <- theta[1]; sigma <- theta[2]; z <- theta[3]
  }

  return (g0 * (1 - exp(-((d/sigma)^-z))))
}

ss.identity <- function (d, theta, cutoff = 130, link = TRUE) {
  b0.ss <- theta[1]
  b1.ss <- theta[2]
  b2.ss <- theta[3]
  sigma.b0.ss <- ifelse(link, exp(theta[4]), theta[4])
  sigma.ss <- ifelse(link, exp(theta[5]), theta[5])

  mean <- b0.ss - b1.ss * d

  return (1 - pnorm(cutoff, mean = mean, sd = sigma.ss))
}

euclideanSq <- function (d, theta.1, theta.2, detfn.1, detfn.2, link = FALSE) {
  (detfn.1(d, theta.1, link = FALSE) - detfn.2(d, theta.2, link = FALSE))^2
}

euclideanNorm <- function(theta.1, theta.2,  detfn.1, detfn.2, link = FALSE) {
  return (pracma::quadinf(euclideanSq, 0, Inf, theta.1 = theta.1, theta.2 = theta.2, detfn.1 = detfn.1, detfn.2 = detfn.2, link = link)$Q^(1/2))
}


# -------------------------------------------------------------------------
palette("Dark 2")
Design.Names <- c("2-by-1 Line", "3-by-1 Line", "4-by-1 Line", "Equil. Triangle", "2-by-2 Grid", "3-by-3 Grid", "5-by-5 Grid")
line.cols = c("hn" = 4, "hr" = 6, "hhn" = 3, "ss" = 1)
line.order = c(1, 2, 4, 3, 5, 6, 7)
Out.Col <- col2rgb(line.cols, alpha = FALSE) |> t() |> rgb(alpha = 100, maxColorValue = 255) |> setNames(names(line.cols))
Box.Col <- line.cols

# fascrNew -------------------------------------------------------------------

fascr.master <- subset(masterList, survey == "fascr")
fascr.master <- fascr.master[!is.na(fascr.master$theta), ]

fascr.D.star <- subset(fascr.master, theta.name == "D") # & theta < 3.000068e+04 )

# fascr.D.star.HR <- subset(fascr.D.star, detfn == "hr")
# fascr.D.star.HR.Info <- aggregate(theta ~ aux, data = fascr.D.star.HR, boxplot.stats, do.conf = FALSE)
#

fascr.res <- expand.grid(Design = unique(fascr.master$design), DetFn = unique(fascr.master$detfn), Aux = unique(fascr.master$aux),
                         relBias = 0, CV = 0, mu = 0, med = 0, sd = 0, stringsAsFactors = FALSE)

for (i in 1:nrow(fascr.res)) {
  inner.x <- subset(fascr.D.star, design == fascr.res$Design[i] & detfn == fascr.res$DetFn[i] & aux == fascr.res$Aux[i])
  fascr.res$relBias[i] <- 100 * (mean(inner.x$theta) - sim.theta$fascr$theta$D) / sim.theta$fascr$theta$D
  fascr.res$mu[i] <- mean(inner.x$theta)
  fascr.res$med[i] <- median(inner.x$theta)
  fascr.res$mse[i] <- mean((inner.x$theta / 25 - sim.theta$fascr$theta$D / 25)^2)
  fascr.res$sd[i] <- sd(inner.x$theta)
  fascr.res$CV[i] <- sd(inner.x$theta) / mean(inner.x$theta) * 100
}

fascr.res.noHR <- subset(fascr.res, DetFn != "hr")
fascr.D.star.noHR <- subset(fascr.D.star, detfn != "hr")

cbind(fascr.res[, c("Design", "DetFn", "Aux")], fascr.res[, c("relBias", "CV", "mse")] |> round(1))


fascr.res[36, 4:9] = NA

# -------------------------------------------------------------------------

fascr.base = subset(masterList, survey == "fascr" & theta.name == "loglik" & detfn != "ss")
fascr.base$design.f = factor(fascr.base$design, levels = c(files)[line.order], labels = Design.Names[line.order])
fascr.seed = split(fascr.base, fascr.base$seed)

boiler.plate = expand.grid(DetFn = names(line.cols), Aux = c(TRUE, FALSE), Design = Design.Names[line.order])
boiler.plate$Freq = 0

toBoil = lapply(fascr.seed, \(x) {
  toCheck = expand.grid(Aux = c(TRUE, FALSE), Design = Design.Names[line.order])
  toCheck$whichDetFn = sapply(1:nrow(toCheck), \(y) {
    z = subset(x, x$aux == toCheck$Aux[y] & as.character(x$design.f) == toCheck$Design[y])
    if (nrow(z) > 0) {
      aics = -2 * z$theta + 2
      for (i in 1:nrow(z)) {
        if (z$detfn[i] == "hn" | z$detfn[i] == "hhn") {
          aics[i] = aics[i] + 2 * 2
        }
        if (z$detfn[i] == "hr") {
          aics[i] = aics[i] + 2 * 3
        }
        if (z$aux[i]) {
          aics[i] = aics[i] + 2 * 1
        }
      }
      names(aics) = z$detfn

      return (which.min(aics) |> names())
    }

    return (NA_character_)
  })

  toCheck}) |>
  do.call(what = "rbind")

toBoil$whichDetFn.f = factor(toBoil$whichDetFn, levels = c("hn", "hr", "hhn"))
toBoil$Design.rev = factor(toBoil$Design, levels = rev(levels(toBoil$Design)),
                           labels = c("5-by-5\nGrid", "3-by-3\nGrid", "2-by-2\nGrid",
                                      "4-by-1\nLine", "Equil.\nTri.", "3-by-1\nLine",
                                      "2-by-1\nLine"))

boiled.vege = xtabs( ~ whichDetFn.f + Design.rev + Aux, data = toBoil)

boiled.vege[3, 7, 1] = 0

palette("Dark 2")

pdf("./figures/05-figureM3.pdf", height = 10.85, width = 11)

layout(matrix(1:2, nrow = 2))

par(mai = c(0.55, 0.45, 0.25, 2.25), cex.axis = 0.8, mgp = c(1.35, 0.35, 0), xaxs = "i", tcl = -0.35)
(100 * prop.table(boiled.vege[, , 2], margin = 2)) |>
  barplot(horiz = TRUE, col = line.cols[1:3], xlab = "AIC Support (%)", xlim = c(0, 100))
mtext("with Times of Arrival", side = 3, font = 2, line = -0.5, cex = 1.25)

par(mai = c(0.55, 0.45, 0.25, 2.25), cex.axis = 0.8, mgp = c(1.35, 0.35, 0), xaxs = "i", tcl = -0.35)
(100 * prop.table(boiled.vege[, , 1], margin = 2)) |>
  barplot(horiz = TRUE, col = line.cols[1:3], xlab = "AIC Support (%)", xlim = c(0, 100))
mtext("without Times of Arrival", side = 3, font = 2, line = -0.5, cex = 1.25)

par(fig = c(0, 1, 0, 0.975), oma = c(0, 0, 0, 0), mar = c(0.25, 0, 0, 0), new = TRUE, cex = 0.9)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bg = "white", legend = c("Halfnormal", "Hazard rate", "Hazard-halfnormal"), pch = 22,
       pt.cex = 2, pt.bg = c(4,6,3), bty = "n", cex = 1.25, title = "Detection function")


dev.off()

# -------------------------------------------------------------------------

palette("Dark 2")


###

pdf("./figures/05-figureM1.pdf", height = 6.5, width = 10)
palette("Dark 2")

layout(matrix(1:2, nrow = 2))
par(mai = c(0.35, 0.5, 0.1, 2.25), cex = 0.8)

# plot(0, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(1, 7), ylim = c(-100, 100))
plot(0, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(1, 7), ylim = range(fascr.res.noHR$relBias))

abline(h = seq(100, -100, by = -20), col = "grey90", lty = 1)

for (i in c("hn", "hhn", "ss")) {
  j = subset(fascr.res, DetFn == i & Aux == FALSE)
  lines(1:7, j$relBias[line.order], type = "b", col = line.cols[i], lty = 6, lwd = 2, cex = 2)

  j = subset(fascr.res, DetFn == i & Aux == TRUE)
  lines(1:7, j$relBias[line.order], type = "b", col = line.cols[i], lty = 1, lwd = 2, cex = 2)
}

axis(2, tck = -0.025, cex.axis = 1.05, mgp = c(1, 0.5, 0)); mtext("Bias (%)", side = 2, line = 1.75)
axis(1, at = 1:7, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0), labels = Design.Names[line.order])
box()

par(mai = c(0.35, 0.5, 0.1, 2.25), cex = 0.8)
plot(0, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(1, 7), ylim = range(fascr.res.noHR$mse))

abline(h = seq(0, 12000, by = 500), col = "grey90", lty = 1)

for (i in c("hn", "hhn", "ss")) {
  j = subset(fascr.res, DetFn == i & Aux == FALSE)
  lines(1:7, j$mse[line.order], type = "b", col = line.cols[i], lty = 6, lwd = 2, cex = 2)

  j = subset(fascr.res, DetFn == i & Aux == TRUE)
  lines(1:7, j$mse[line.order], type = "b", col = line.cols[i], lty = 1, lwd = 2, cex = 2)
}

axis(2, tck = -0.025, cex.axis = 1.05, mgp = c(1, 0.5, 0)); mtext("MSE", side = 2, line = 1.75)
axis(1, at = 1:7, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0), labels = Design.Names[line.order])
box()

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, xaxs = "i", yaxs = "i", cex = 0.9)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', ylim = c(0 , 1), xlim = c(0 , 1))
legend(x = 0.79, y = 0.675, legend = c("Halfnormal", "Hazard-halfnormal", "Signal Strength"),
       col = line.cols[c(1, 3, 4)], lwd = 2, cex = 1, title = "Detection Function", bty = "n")

legend(x = 0.79, y = 0.49, lty = c(1, 6), legend = c("Included Submodel", "Excluded Submodel"),
       lwd = 2, cex = 1, title = "Times of Arrival", bty = "n")

dev.off()

##
palette("Dark 2")

fascr.D.star$design.f = factor(fascr.D.star$design, levels = c(files)[line.order])
fascr.D.star$theta[which(fascr.D.star$detfn == "hhn" & fascr.D.star$aux == FALSE & fascr.D.star$design == "fascrNew_1.csv")] = NA
### Mistake picked up here...

# pdf("./figures/05-figureM2.pdf", height = 12, width = 8)
pdf("./figures/05-figureM2.pdf", height = 8, width = 16)

layout(matrix(c(1,3,2,4), nrow = 2, ncol = 2))
# layout(matrix(1:4, nrow = 4))
par(mai = c(0.3, 0.5, 0.1, 0.25), cex = 0.8)

## Halfnormal
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "hn" & aux == TRUE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["hn"], col = Box.Col["hn"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 - 0.25, xlim = c(0.5, 7.5), whisklty = 1, staplewex = 0.75, medlwd = 2)
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "hn" & aux == FALSE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["hn"], col = Out.Col["hn"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 + 0.25, add = TRUE, whisklty = 1, staplewex = 0.75, medlwd = 2)
legend("topright", legend = c("w/ Times of Arrival", "w/out Times of Arrival"), pch = 22, bty = "n", pt.cex = 2,
       pt.bg = c(Box.Col["hn"], Out.Col["hn"]))

axis(2, tck = -0.025, cex.axis = 1.05, mgp = c(1, 0.5, 0))
title(ylab = expression("Call density p/ha p/s"), las = 1, line = 1.5, cex.lab = 1.25)
axis(1, at = 1:7, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0), labels = Design.Names[line.order])

abline(h = sim.theta$fascr$theta$D / 25, col = "black", lwd = 3)
abline(h = sim.theta$fascr$theta$D / 25, col = "red", lwd = 1)

with(subset(fascr.res, DetFn == "hn" & Aux == TRUE), {points(x = line.order - 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})
with(subset(fascr.res, DetFn == "hn" & Aux == FALSE), {points(x = line.order + 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})

text(x = 7.9, y = diff(range(fascr.D.star.noHR$theta)) / 25 / 2, xpd = TRUE, labels = "Halfnormal", font = 2, srt = 270, cex = 1.5)
box()

## Hazard rate
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "hr" & aux == TRUE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["hr"], col = Box.Col["hr"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 - 0.25, xlim = c(0.5, 7.5), whisklty = 1, staplewex = 0.75, medlwd = 2)
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "hr" & aux == FALSE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["hr"], col = Out.Col["hr"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 + 0.25, add = TRUE, whisklty = 1, staplewex = 0.75, medlwd = 2)
legend("topright", legend = c("w/ Times of Arrival", "w/out Times of Arrival"), pch = 22, bty = "n", pt.cex = 2,
       pt.bg = c(Box.Col["hr"], Out.Col["hr"]))

axis(2, tck = -0.025, cex.axis = 1.05, mgp = c(1, 0.5, 0))
title(ylab = expression("Call density p/ha p/s"), las = 1, line = 1.5, cex.lab = 1.25)
axis(1, at = 1:7, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0), labels = Design.Names[line.order])

abline(h = sim.theta$fascr$theta$D / 25, col = "black", lwd = 3)
abline(h = sim.theta$fascr$theta$D / 25, col = "red", lwd = 1)

with(subset(fascr.res, DetFn == "hr" & Aux == TRUE), {points(x = line.order - 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})
with(subset(fascr.res, DetFn == "hr" & Aux == FALSE), {points(x = line.order + 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})

text(x = 7.9, y = diff(range(fascr.D.star.noHR$theta)) / 25 / 2, xpd = TRUE, labels = "Hazard rate", font = 2, srt = 270, cex = 1.5)
box()

## Hazard halfnormal
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "hhn" & aux == TRUE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["hhn"], col = Box.Col["hhn"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 - 0.25, xlim = c(0.5, 7.5), whisklty = 1, staplewex = 0.75, medlwd = 2)
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "hhn" & aux == FALSE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["hhn"], col = Out.Col["hhn"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 + 0.25, add = TRUE, whisklty = 1, staplewex = 0.75, medlwd = 2)
legend("topright", legend = c("w/ Times of Arrival", "w/out Times of Arrival"), pch = 22, bty = "n", pt.cex = 2,
       pt.bg = c(Box.Col["hhn"], Out.Col["hhn"]))

abline(h = sim.theta$fascr$theta$D / 25, col = "black", lwd = 3)
abline(h = sim.theta$fascr$theta$D / 25, col = "red", lwd = 1)


axis(2, tck = -0.025, cex.axis = 1.05, mgp = c(1, 0.5, 0))
title(ylab = expression("Call density p/ha p/s"), las = 1, line = 1.5, cex.lab = 1.25)
axis(1, at = 1:7, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0), labels = Design.Names[line.order])

with(subset(fascr.res, DetFn == "hhn" & Aux == TRUE), {points(x = line.order - 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})
with(subset(fascr.res, DetFn == "hhn" & Aux == FALSE), {points(x = line.order + 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})

text(x = 7.9, y = diff(range(fascr.D.star.noHR$theta)) / 25 / 2, xpd = TRUE, labels = "Hazard-halfnormal", font = 2, srt = 270, cex = 1.5)
box()

## Signal strength
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "ss" & aux == TRUE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["ss"], col = Box.Col["ss"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 - 0.25, xlim = c(0.5, 7.5), whisklty = 1, staplewex = 0.75, medlwd = 2)
boxplot(theta / 25 ~ design.f, data = subset(fascr.D.star, detfn == "ss" & aux == FALSE), boxwex = 0.3, xlab = "", ylab = "", axes = FALSE,
        outpch = 16, outcol = Out.Col["ss"], col = Out.Col["ss"], ylim = range(fascr.D.star.noHR$theta) / 25,
        at = 1:7 + 0.25, add = TRUE, whisklty = 1, staplewex = 0.75, medlwd = 2)
legend("topright", legend = c("w/ Times of Arrival", "w/out Times of Arrival"), pch = 22, bty = "n", pt.cex = 2,
       pt.bg = c(Box.Col["ss"], Out.Col["ss"]))

axis(2, tck = -0.025, cex.axis = 1.05, mgp = c(1, 0.5, 0));
title(ylab = expression("Call density p/ha p/s"), las = 1, line = 1.5, cex.lab = 1.25)
axis(1, at = 1:7, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0), labels = Design.Names[line.order])

abline(h = sim.theta$fascr$theta$D / 25, col = "black", lwd = 3)
abline(h = sim.theta$fascr$theta$D / 25, col = "red", lwd = 1)

with(subset(fascr.res, DetFn == "ss" & Aux == TRUE), {points(x = line.order - 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})
with(subset(fascr.res, DetFn == "ss" & Aux == FALSE), {points(x = line.order + 0.25, y = mu / 25, bg = "white", pch = 23, cex = 1.25)})

text(x = 7.9, y = diff(range(fascr.D.star.noHR$theta)) / 25 / 2, xpd = TRUE, labels = "Signal strength", font = 2, srt = 270, cex = 1.5)
box()

dev.off()

##
# Hazard rate untruncated
