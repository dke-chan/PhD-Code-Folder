load("./data/sonobuoys.RData")
source("./code/sonobuoys/00-aux-funcs.R")
library(ascr); library(xtable); library(tidyr); library(dplyr)

ascr.models <- readRDS("./code/sonobuoys/ascr-models.RDS")
logistic.like.models <- readRDS("./code/sonobuoys/logistic-like-models.RDS")
mrds.models <- readRDS("./code/sonobuoys/mrds-models.RDS")
logistic.like.models$spec <- "Logistic"
mrds.models$spec <- "MRDS"

# Standard is unchanged (otherwise, 1 / 24)
# Batched is unchanged (otherwise, 10 / 24)
# Unbatched is density is multiplied by 23 / 24

# -------------------------------------------------------------------------

## (Line 73)
surveyLength.session <- seaTrial.df %>%
  select(TrackedGroup, posixDate) %>%
  mutate(TrackedGroup = LETTERS[as.numeric(gsub("S([0-9]+)", "\\1", TrackedGroup))]) %>%
  filter(complete.cases(.)) %>%
  mutate(dateTime = lubridate::parse_date_time(posixDate, "dmYHM"),
         minute = lubridate::minute(dateTime)) %>%
  group_by(TrackedGroup) %>%
  summarise(duration = max(minute) - min(minute)) %>%
  mutate(duration.new = ifelse(duration == 0, 1, duration))

# Table 3.1 ---------------------------------------------------------------
ascr.toPrint <- ascr.models %>%
  mutate(spec.fact = factor(spec, levels = c("start", "start.abs", "drift.batch", "drift.unbatched")),
         detfn.fact = factor(detfn, levels = c("hn", "hr", "th", "hhn")),
         est.mod = if_else(est.name == "D", est * 100 * 100, est),
         # est.mod = if_else(est.name == "D" & spec.fact == "start", est.mod / 24, est.mod),
         toReport.A = sprintf("%.1f", est.mod),
         se = if_else(est.name == "D", se * 100 * 100, est),
         # se = if_else(est.name == "D"  & spec.fact == "start", se / 24, se),
         toReport.B = paste0("(", sprintf("%.2f", se), ")")) %>%
  filter(spec.fact != "start.abs") %>%
  arrange(detfn.fact, spec.fact) %>%
  select(detfn, est.name, toReport.A, toReport.B, spec) %>%
  pivot_wider(id_cols = c(detfn, est.name), values_from = c(toReport.A, toReport.B), names_from = c(spec)) %>%
  as.data.frame()

ascr.toPrint.order <- sapply(c("start", "start.abs", "drift.batch", "drift.unbatched"), \(x) {
  paste0(c("toReport.A_", "toReport.B_"), x)
}) |> c()

ascr.toPrint <- ascr.toPrint[, c("detfn", "est.name", ascr.toPrint.order)]

ascr.toPrint$est.name <- sapply(ascr.toPrint$est.name, \(x) {
  switch(x, D = "$\\widehat{D}$", g0 = "$\\widehat{g_0}$", sigma = "$\\widehat{\\sigma}$", lambda0 = "$\\widehat{\\lambda_0}$",
         z = "$\\widehat{z}$", shape = "$\\widehat{\\alpha}$", scale = "$\\widehat{\\lambda}$",
         kappa = "$\\widehat{\\kappa}$")
})

ascr.toPrint$detfn <- sapply(ascr.toPrint$detfn, \(x) {
  switch(x, hn = "Halfnormal", hr = "Hazard rate", th = "Threshold", hhn = "Hazard halfnormal")
})
ascr.toPrint$detfn[(1:17)[-c(1, 5, 10, 14)]] <- ""

colnames(ascr.toPrint)[1:2] <- c("$g\\{\\bm{s}, \\bm{x}_k(\\tau_i); \\bm{\\theta}\\}$", "")

ascr.Table <- xtable(ascr.toPrint, caption = "toInsert", label = "toInsert")

print(ascr.Table, sanitize.text.function = \(x){x}, booktabs = TRUE, caption.placement = "bottom",
      tabular.environment = "tabular", include.rownames = FALSE)

# Table 3.2 ---------------------------------------------------------------
combined.toPrint <- bind_rows(mrds.models, logistic.like.models) %>%
  mutate(spec.fact = factor(spec, levels = c("MRDS", "Logistic")),
         detfn.fact = factor(detfn, levels = c("hn", "hr", "th", "hhn")),
         est.mod = if_else(est.name == "D", est * 100 * 100, est),
         se = if_else(est.name == "D", se * 100 * 100, est),
         toReport.A = sprintf("%.2f", est.mod),
         toReport.B = paste0("(", sprintf("%.2f", se), ")")) %>%
  arrange(detfn.fact, spec.fact) %>%
  select(detfn, est.name, toReport.A, toReport.B, spec) %>%
  pivot_wider(id_cols = c(detfn, est.name), values_from = c(toReport.A, toReport.B), names_from = c(spec)) %>%
  as.data.frame()

combined.toPrint.order <- sapply(c("MRDS", "Logistic"), \(x) {
  paste0(c("toReport.A_", "toReport.B_"), x)
}) |> c()

combined.toPrint <- combined.toPrint[, c("detfn", "est.name", combined.toPrint.order)]

combined.toPrint$est.name <- sapply(combined.toPrint$est.name, \(x) {
  switch(x, D = "$\\widehat{D}$", g0 = "$\\widehat{g_0}$", sigma = "$\\widehat{\\sigma}$", lambda0 = "$\\widehat{\\lambda_0}$",
         z = "$\\widehat{z}$", shape = "$\\widehat{\\alpha}$", scale = "$\\widehat{\\lambda}$",
         kappa = "$\\widehat{\\kappa}$")
})

combined.toPrint$detfn <- sapply(combined.toPrint$detfn, \(x) {
  switch(x, hn = "Halfnormal", hr = "Hazard rate", th = "Threshold", hhn = "Hazard halfnormal")
})
combined.toPrint$detfn[(1:17)[-c(1, 5, 10, 14)]] <- ""

colnames(combined.toPrint)[1:2] <- c("$g\\{\\bm{s}, \\bm{x}_k(\\tau_i); \\bm{\\theta}\\}$", "")

combined.Table <- xtable(combined.toPrint, caption = "toInsert", label = "toInsert")

print(combined.Table, sanitize.text.function = \(x){x}, booktabs = TRUE, caption.placement = "bottom",
      tabular.environment = "tabular", include.rownames = FALSE, NA.string = "--")

# Figure 3.6 --------------------------------------------------------------

layout.mat <- c(rep(rep(c(1, 2), each = 8), times = 8),
  rep(rep(c(3, 4), each = 8), times = 8)) |>
  matrix(nrow = 16, ncol = 16, byrow = TRUE) |>
  rbind(rep(5, 16))

model.cols <- RColorBrewer::brewer.pal(3, "Dark2")
d <- 0:5250

pdf("./figures/03-figure06.pdf", height = 9, width = 8)
par(cex = 1)
layout(layout.mat)
par(mai = c(0.45, 0.5, 0.3, 0.05), xaxs = "i")

## HN
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Halfnormal", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = halfnormal(d, subset(ascr.models, detfn == "hn" & spec == "drift.unbatched")$est[2:3], link = FALSE),
      col = model.cols[1], lwd = 2)
lines(x = d, y = halfnormal(d, subset(mrds.models, detfn == "hn" & spec == "MRDS")$est[2:3], link = FALSE),
      col = model.cols[2], lwd = 2)
lines(x = d, y = halfnormal(d, subset(logistic.like.models, detfn == "hn" & spec == "Logistic")$est, link = FALSE),
      col = model.cols[3], lwd = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## HR
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Hazard rate", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = hazardRate(d, subset(ascr.models, detfn == "hr" & spec == "drift.unbatched")$est[2:4], link = FALSE),
      col = model.cols[1], lwd = 2)
lines(x = d, y = hazardRate(d, subset(mrds.models, detfn == "hr" & spec == "MRDS")$est[2:4], link = FALSE),
      col = model.cols[2], lwd = 2)
lines(x = d, y = hazardRate(d, subset(logistic.like.models, detfn == "hr" & spec == "Logistic")$est, link = FALSE),
      col = model.cols[3], lwd = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## TH
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Threshold", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = threshold(d, subset(ascr.models, detfn == "th" & spec == "drift.unbatched")$est[2:3], link = FALSE),
      col = model.cols[1], lwd = 2)
lines(x = d, y = threshold(d, subset(mrds.models, detfn == "th" & spec == "MRDS")$est[2:3], link = FALSE),
      col = model.cols[2], lwd = 2)
lines(x = d, y = threshold(d, subset(logistic.like.models, detfn == "th" & spec == "Logistic")$est, link = FALSE),
      col = model.cols[3], lwd = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## HHN
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Hazard-halfnormal", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = hazardHalfnormal(d, subset(ascr.models, detfn == "hhn" & spec == "drift.unbatched")$est[2:3], link = FALSE),
      col = model.cols[1], lwd = 2)
lines(x = d, y = hazardHalfnormal(d, subset(mrds.models, detfn == "hhn" & spec == "MRDS")$est[2:3], link = FALSE),
      col = model.cols[2], lwd = 2)
lines(x = d, y = hazardHalfnormal(d, subset(logistic.like.models, detfn == "hhn" & spec == "Logistic")$est, link = FALSE),
      col = model.cols[3], lwd = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## Legend
par(mai = c(0.05, 0.05, 0.05, 0.05), yaxs = "r")
plot(0, xlim = c(0, 10), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")

# text(5, 0.6, labels = "Legend", pos = 3, cex = 1.25, font = 2)

segments(2.5, 0.5, 3, lwd = 2, col = model.cols[1])
text(3, 0.475, labels = "Acoustic SCR", pos = 4, cex = 1.25)

segments(4.5, 0.5, 5, lwd = 2, col = model.cols[2])
text(5, 0.475, labels = "Acoustic MRDS", pos = 4, cex = 1.25)

segments(6.5, 0.5, 7, lwd = 2, col = model.cols[3])
text(7, 0.475, labels = "Bernoulli", pos = 4, cex = 1.25)
dev.off()

# Table 3.3 ---------------------------------------------------------------

obs.euclideans <- lapply(c("hn", "hr", "th", "hhn"), \(x) {
  detfn.param <- switch(x, hn = c("g0", "sigma"), hr = c("g0", "sigma", "z"),
                        th = c("shape", "scale"), hhn = c("lambda0", "sigma"))
  detfn.f <- switch(x, hn = halfnormal, hr = hazardRate,
                    th = threshold, hhn = hazardHalfnormal)

  a <- subset(ascr.models, detfn == x)
  m <- subset(mrds.models, detfn == x)
  l <- subset(logistic.like.models, detfn == x)

  L2Norms <- c(euclideanNorm(a$est[a$est.name %in% detfn.param], l$est, detfn.f, link = FALSE),
               euclideanNorm(m$est[m$est.name %in% detfn.param], l$est, detfn.f, link = FALSE),
               euclideanNorm(a$est[a$est.name %in% detfn.param], m$est[m$est.name %in% detfn.param], detfn.f, link = FALSE))
}) |> do.call(what = "rbind") |> as.data.frame()

colnames(obs.euclideans) <- c("A-L", "M-L", "A-M")

obs.euclideans$detfn <- c("hn", "hr", "th", "hhn")

xtable(obs.euclideans[, c(4, 1:3)], caption = "toInsert", label = "toInsert") |>
  print(caption.placement = "bottom", tabular.environment = "tabular", include.rownames = FALSE, booktabs = TRUE)

# Figure 3.7 --------------------------------------------------------------
model.cols.2 = c("#e41a1c", "#377eb8", "#4daf4a")

transparent.cols <- model.cols.2 |> col2rgb() |> t() |> rgb(alpha = 175, maxColorValue = 255)
hn.SimStudy <- read.csv("./code/sonobuoys/hnSim.csv")
hr.SimStudy <- read.csv("./code/sonobuoys/hrSim.csv")
th.SimStudy <- read.csv("./code/sonobuoys/thSim.csv")
hhn.SimStudy <- read.csv("./code/sonobuoys/hhnSim.csv")

# layout(matrix(1:8, nrow = 4, ncol = 2, byrow = 2))

hn.sim.L2Norms <- subset(hn.SimStudy, whichOutput == "A-L" | whichOutput == "M-L" | whichOutput == "A-M")
hn.sim.L2Norms$whichFactor <- factor(hn.sim.L2Norms$whichOutput, levels = c("A-M", "A-L", "M-L"), labels = c("A-M hn", "A-L hn", "M-L hn"))
hr.sim.L2Norms <- subset(hr.SimStudy, whichOutput == "A-L" | whichOutput == "M-L" | whichOutput == "A-M")
hr.sim.L2Norms$whichFactor <- factor(hr.sim.L2Norms$whichOutput, levels = c("A-M", "A-L", "M-L"), labels = c("A-M hr", "A-L hr", "M-L hr"))
th.sim.L2Norms <- subset(th.SimStudy, whichOutput == "A-L" | whichOutput == "M-L" | whichOutput == "A-M")
th.sim.L2Norms$whichFactor <- factor(th.sim.L2Norms$whichOutput, levels = c("A-M", "A-L", "M-L"), labels = c("A-M th", "A-L th", "M-L th"))
hhn.sim.L2Norms <- subset(hhn.SimStudy, whichOutput == "A-L" | whichOutput == "M-L" | whichOutput == "A-M")
hhn.sim.L2Norms$whichFactor <- factor(hhn.sim.L2Norms$whichOutput, levels = c("A-M", "A-L", "M-L"), labels = c("A-M hhn", "A-L hhn", "M-L hhn"))

pdf("./figures/03-figure07.pdf", height = 5.5, width = 12)
par(mai = c(0.6, 0.475, 0.1, 0.1), xaxs = "i")

boxplot(output ~ whichFactor, data = rbind(hn.sim.L2Norms, hr.sim.L2Norms, th.sim.L2Norms, hhn.sim.L2Norms),
        boxwex = 0.3, col = model.cols.2[1:3], outpch = 16, outcol = transparent.cols[1:3], axes = FALSE, xlab = "", ylab = "")
abline(v = seq(3.5, 9.5, by = 3), col = "grey75")

points(x = seq(1, 10, by = 3), y = obs.euclideans$`A-M`, pch = 23, cex = 1.5, bg = "ghostwhite", lwd = 2)
points(x = seq(2, 11, by = 3), y = obs.euclideans$`A-L`, pch = 23, cex = 1.5, bg = "ghostwhite", lwd = 2)
points(x = seq(3, 12, by = 3), y = obs.euclideans$`M-L`, pch = 23, cex = 1.5, bg = "ghostwhite", lwd = 2)
# points(x = seq(2, 8, by = 2), y = obs.euclideans$`A-M`, pch = 23, cex = 1, bg = "ghostwhite")

axis(1, at = 1:12, labels = rep(c("SCR and\nMRDS", "SCR and\nBernoulli", "MRDS and\nBernoulli"), 4), tck = -0.015, mgp = c(1, 1.5, 0), cex.axis = 1)
# axis(1, at = 1:8, labels = rep(c("SCR - Bernoulli", "MRDS - Bernoulli"), 4), tck = -0.015, mgp = c(1, 0.2, 0), cex.axis = 0.7)
axis(3, at = seq(2, 11, by = 3), labels = c("Halfnormal", "Hazard rate", "Threshold", "Hazard-halfnormal"), tck = 0, cex.axis = 1.1, font = 2, mgp = c(-1, -1.5, 0))
axis(2, tck = -0.015, mgp = c(1, 0.5, 0), cex.axis = 0.8); mtext("Euclidean norm", side = 2, line = 1.35, font = 1, cex.axis = 1)
box()
dev.off()


# Table 3.3 Extra ---------------------------------------------------------

# mean(subset(hn.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[1])
# mean(subset(hn.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[1])
# mean(subset(hr.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[2])
# mean(subset(hr.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[2])
# mean(subset(th.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[3])
# mean(subset(th.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[3])
# mean(subset(hhn.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[4])
# mean(subset(hhn.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[4])

mean(subset(hn.sim.L2Norms, whichOutput == "A-M")$output > obs.euclideans$`A-M`[1])
mean(subset(hn.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[1])
mean(subset(hn.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[1])
mean(subset(hr.sim.L2Norms, whichOutput == "A-M")$output > obs.euclideans$`A-M`[2])
mean(subset(hr.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[2])
mean(subset(hr.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[2])
mean(subset(th.sim.L2Norms, whichOutput == "A-M")$output > obs.euclideans$`A-M`[3])
mean(subset(th.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[3])
mean(subset(th.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[3])
mean(subset(hhn.sim.L2Norms, whichOutput == "A-M")$output > obs.euclideans$`A-M`[4])
mean(subset(hhn.sim.L2Norms, whichOutput == "A-L")$output > obs.euclideans$`A-L`[4])
mean(subset(hhn.sim.L2Norms, whichOutput == "M-L")$output > obs.euclideans$`M-L`[4])

# Figure 3.8 --------------------------------------------------------------
hn.sim.ndet <- subset(hn.SimStudy, whichOutput == "n.detect")$output
hr.sim.ndet <- subset(hr.SimStudy, whichOutput == "n.detect")$output
th.sim.ndet <- subset(th.SimStudy, whichOutput == "n.detect")$output
hhn.sim.ndet <- subset(hhn.SimStudy, whichOutput == "n.detect")$output

x <- halfnormal(full.distance.mat * 1000, logistic.like.models$est[which(logistic.like.models$detfn == "hn")], link = FALSE)
hn.est <- sum(1 - apply(1 - x, 1, prod))
x <- hazardRate(full.distance.mat * 1000, logistic.like.models$est[which(logistic.like.models$detfn == "hr")], link = FALSE)
hr.est <- sum(1 - apply(1 - x, 1, prod))
x <- threshold(full.distance.mat * 1000, logistic.like.models$est[which(logistic.like.models$detfn == "th")], link = FALSE)
th.est <- sum(1 - apply(1 - x, 1, prod))
x <- hazardHalfnormal(full.distance.mat * 1000, logistic.like.models$est[which(logistic.like.models$detfn == "hhn")], link = FALSE)
hhn.est <- sum(1 - apply(1 - x, 1, prod))

pdf("./figures/03-figure08.pdf", height = 5, width = 9)
par(mai = c(0.4, 0.475, 0.1, 0.1))
boxplot(cbind(hn.sim.ndet, hr.sim.ndet, th.sim.ndet, hhn.sim.ndet), boxwex = 0.25, col = "lightblue",
        outpch = 21, outbg = "lightblue", xlab = "", ylab = "", axes = FALSE)
abline(h = 23, lty = 5, col = "red")
boxplot(cbind(hn.sim.ndet, hr.sim.ndet, th.sim.ndet, hhn.sim.ndet), boxwex = 0.25, col = "lightblue",
        outpch = 21, outbg = "lightblue", xlab = "", ylab = "", axes = FALSE, add = TRUE)

points(1:4, y = c(hn.est, hr.est, th.est, hhn.est), pch = 23, cex = 1.25, bg = "ghostwhite")

axis(1, tck = -0.015, at= 1:4, labels = c("Halfnormal", "Hazard rate", "Threshold", "Hazard-halfnormal"), mgp = c(1, 0.5, 0), cex.axis = 1, font = 2)#; mtext("Detection function", side = 1, line = 1.325, font = 2, cex = 1)
axis(2, tck = -0.015, mgp = c(1, 0.5, 0), cex.axis = 0.8); mtext("Number of detections", side = 2, line = 1.35, font = 2, cex = 1.1)

box()
dev.off()

# -------------------------------------------------------------------------

hn.sim.failedit <- subset(hn.SimStudy, whichOutput == "failed.it")$output
hr.sim.failedit <- subset(hr.SimStudy, whichOutput == "failed.it")$output
th.sim.failedit <- subset(th.SimStudy, whichOutput == "failed.it")$output
hhn.sim.failedit <- subset(hhn.SimStudy, whichOutput == "failed.it")$output


# Figure Bonus ------------------------------------------------------------

palette("Dark 2")

pdf("./figures/03-figureBB.pdf", height = 9, width = 8)
par(cex = 1)
layout(layout.mat)
par(mai = c(0.45, 0.5, 0.3, 0.05), xaxs = "i")

## HN
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Halfnormal", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = halfnormal(d, subset(ascr.models, detfn == "hn" & spec == "start")$est[2:3], link = FALSE),
      col = 4, lwd = 2, lty = 1)
lines(x = d, y = halfnormal(d, subset(ascr.models, detfn == "hn" & spec == "drift.batch")$est[2:3], link = FALSE),
      col = 6, lwd = 2, lty = 5)
lines(x = d, y = halfnormal(d, subset(ascr.models, detfn == "hn" & spec == "drift.unbatched")$est[2:3], link = FALSE),
      col = 1, lwd = 2, lty = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## HR
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Hazard rate", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = hazardRate(d, subset(ascr.models, detfn == "hr" & spec == "start")$est[2:4], link = FALSE),
      col = 4, lwd = 2, lty = 1)
lines(x = d, y = hazardRate(d, subset(ascr.models, detfn == "hr" & spec == "drift.batch")$est[2:4], link = FALSE),
      col = 6, lwd = 2, lty = 5)
lines(x = d, y = hazardRate(d, subset(ascr.models, detfn == "hr" & spec == "drift.unbatched")$est[2:4], link = FALSE),
      col = 1, lwd = 2, lty = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## TH
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Threshold", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = threshold(d, subset(ascr.models, detfn == "th" & spec == "start")$est[2:3], link = FALSE),
      col = 4, lwd = 2, lty = 1)
lines(x = d, y = threshold(d, subset(ascr.models, detfn == "th" & spec == "drift.batch")$est[2:3], link = FALSE),
      col = 6, lwd = 2, lty = 5)
lines(x = d, y = threshold(d, subset(ascr.models, detfn == "th" & spec == "drift.unbatched")$est[2:3], link = FALSE),
      col = 1, lwd = 2, lty = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## HHN
plot(0, xlim = c(0, 5200), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
mtext("Hazard-halfnormal", side = 3, line = 0.25, font = 2)
abline(h = c(0, 1), col = "grey75")

lines(x = d, y = hazardHalfnormal(d, subset(ascr.models, detfn == "hhn" & spec == "start")$est[2:3], link = FALSE),
      col = 4, lwd = 2, lty = 1)
lines(x = d, y = hazardHalfnormal(d, subset(ascr.models, detfn == "hhn" & spec == "drift.batch")$est[2:3], link = FALSE),
      col = 6, lwd = 2, lty = 5)
lines(x = d, y = hazardHalfnormal(d, subset(ascr.models, detfn == "hhn" & spec == "drift.unbatched")$est[2:3], link = FALSE),
      col = 1, lwd = 2, lty = 2)

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.55, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()

## Legend
par(mai = c(0.05, 0.05, 0.05, 0.05), yaxs = "r")
plot(0, xlim = c(0, 10), ylim = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")

# text(5, 0.6, labels = "Legend", pos = 3, cex = 1.25, font = 2)

segments(2.5, 0.5, 3, lwd = 2, col = 4, lty = 1)
text(3, 0.475, labels = "Start (Relative)", pos = 4, cex = 1.25)

segments(4.5, 0.5, 5, lwd = 2, col = 6, lty = 5)
text(5, 0.475, labels = "Batched", pos = 4, cex = 1.25)

segments(6.5, 0.5, 7, lwd = 2, col = 1, lty = 2)
text(7, 0.475, labels = "Unbatched", pos = 4, cex = 1.25)
dev.off()


# Figure 3.9 --------------------------------------------------------------

# start.hhn <- fit.ascr(capt = collapsed.capt, traps = buoy.start.xy,
#                       mask = mask.start, detfn = "hhn")
# start.abs.hhn <- fit.ascr(capt = collapsed.capt, traps = buoy.abs.start.xy,
#                           mask = mask.abs.start, detfn = "hhn")
#
# pdf("./figures/03-figure09.pdf", height = 9, width = 8)
# layout(matrix(1:4, 2, 2))
# par(mai = c(0.55, 0.55, 0.3, 0.1))
#
# #
# plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
#      xlim = range(c(seaTrial.df$buoy.utm.x, seaTrial.df$boat.utm.x)) + 500 * c(-1, 1),
#      ylim = range(c(seaTrial.df$buoy.utm.y, seaTrial.df$boat.utm.y)) + 500 * c(-1, 1))
#
# locations(start.abs.hn, 1, add = TRUE, show.labels = FALSE, trap.col = "tomato",
#           cols = list(combined = "steelblue"), lty = list(combined = "longdash"))
#
# points(buoy.abs.start.xy, pch = 21, bg = "tomato", cex = 1)
# seaTrial.df %>% filter(PlaybackNumber == 1) %>% dplyr::select(boat.utm.x, boat.utm.y) %>%
#   points(pch = 21, bg = "steelblue", cex = 1.5)
#
# axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 11S)", side = 1, line = 1.55, font = 2, cex = 0.75)
# axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 11S)", side = 2, line = 1.55, font = 2, cex = 0.8)
#
# # legend("topright", legend = c("First detected call"), cex = 2, bty = "n")
# mtext("First detected call", side = 3, line = 0.25, font = 2)
#
# box()
#
# plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
#      xlim = range(c(seaTrial.df$buoy.utm.x, seaTrial.df$boat.utm.x)) + 500 * c(-1, 1),
#      ylim = range(c(seaTrial.df$buoy.utm.y, seaTrial.df$boat.utm.y)) + 500 * c(-1, 1))
#
# locations(start.abs.hn, 23, add = TRUE, show.labels = FALSE, trap.col = "tomato",
#           cols = list(combined = "steelblue"), lty = list(combined = "longdash"))
#
# points(buoy.abs.start.xy, pch = 21, bg = "tomato", cex = 1)
# seaTrial.df %>% filter(PlaybackNumber == 26) %>% dplyr::select(boat.utm.x, boat.utm.y) %>%
#   points(pch = 21, bg = "steelblue", cex = 1.5)
#
# axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (UTM Zone 11S)", side = 1, line = 1.55, font = 2, cex = 0.75)
# axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (UTM Zone 11S)", side = 2, line = 1.55, font = 2, cex = 0.8)
#
# # legend("topright", legend = c("Last detected call"), cex = 2, bty = "n")
# mtext("Last detected call", side = 3, line = 0.25, font = 2)
#
# box()
#
# #
# plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
#      xlim = range(c(seaTrial.df$buoy.rel.x, seaTrial.df$boat.rel.x)) + 500 * c(-1, 1),
#      ylim = range(c(seaTrial.df$buoy.rel.y, seaTrial.df$boat.rel.y)) + 500 * c(-1, 1))
#
# locations(start.hn, 1, add = TRUE, show.labels = FALSE, trap.col = "tomato",
#           cols = list(combined = "steelblue"), lty = list(combined = "longdash"))
#
# points(buoy.start.xy, pch = 21, bg = "tomato", cex = 1)
# seaTrial.df %>% filter(PlaybackNumber == 1) %>% dplyr::select(boat.rel.x, boat.rel.y) %>%
#   points(pch = 21, bg = "steelblue", cex = 1.5)
#
# axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (Relative)", side = 1, line = 1.55, font = 2, cex = 0.75)
# axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (Relative)", side = 2, line = 1.55, font = 2, cex = 0.8)
#
# # legend("topright", legend = c("First detected call"), cex = 2, bty = "n")
# mtext("First detected call", side = 3, line = 0.25, font = 2)
#
# box()
#
# plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1,
#      xlim = range(c(seaTrial.df$buoy.rel.x, seaTrial.df$boat.rel.x)) + 500 * c(-1, 1),
#      ylim = range(c(seaTrial.df$buoy.rel.y, seaTrial.df$boat.rel.y)) + 500 * c(-1, 1))
#
# locations(start.hn, 23, add = TRUE, show.labels = FALSE, trap.col = "tomato",
#           cols = list(combined = "steelblue"), lty = list(combined = "longdash"))
#
# points(buoy.start.xy, pch = 21, bg = "tomato", cex = 1)
# seaTrial.df %>% filter(PlaybackNumber == 26) %>% dplyr::select(boat.rel.x, boat.rel.y) %>%
#   points(pch = 21, bg = "steelblue", cex = 1.5)
#
# axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("X (Relative)", side = 1, line = 1.55, font = 2, cex = 0.75)
# axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Y (Relative)", side = 2, line = 1.55, font = 2, cex = 0.8)
#
# # legend("topright", legend = c("Last detected call"), cex = 2, bty = "n")
# mtext("Last detected call", side = 3, line = 0.25, font = 2)
#
# box()
#
# dev.off()
