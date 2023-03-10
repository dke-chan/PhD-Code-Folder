library(ascr)

# reference ---------------------------------------------------------------
# identity
## mu_ss = detpars(1) - detpars(2)*column(dists(s), i);
# log
## mu_ss = exp(detpars(1) - detpars(2)*column(dists(s), i));
# spherical
## mu_ss = detpars(1) - 10*log10(square(column(dists(s), i))) - detpars(2)*(column(dists(s), i) - 1.0);

## out <- 1 - pnorm(cutoff, mean = mean, sd = sigma.ss)

# frogs (admbsecr) --------------------------------------------------------
load(file = "code/miss-v2/lightfooti.RData")
lightfooti.admb = lightfooti
fcapt <- lightfooti.admb$capt.trunc
fmask <- lightfooti.admb$mask
ftraps <- lightfooti.admb$traps

# Default settings get the estimates in Stevenson et al. (2015) model estimates in Table 1
fascr <- fit.ascr(fcapt, ftraps, fmask, detfn = "ss",
                  ss.opts = list(cutoff = 140, ss.link = "identity"))
## The density estimate is divided by 25 in Stevenson et al. (2015) to convert it
## density per hectare per second

fascr.log <- fit.ascr(fcapt, ftraps, fmask, detfn = "ss",
                      ss.opts = list(cutoff = 140, ss.link = "log"))
fascr.sphere <- fit.ascr(fcapt, ftraps, fmask, detfn = "ss",
                         ss.opts = list(cutoff = 140, ss.link = "spherical"))

# frog hn and hhn models --------------------------------------------------
fcapt.nullSS <- list(bincapt = fcapt$bincapt, toa = fcapt$toa)
fascr.hn <- fit.ascr(fcapt.nullSS, ftraps, fmask, detfn = "hn")
fascr.hhn <- fit.ascr(fcapt.nullSS, ftraps, fmask, detfn = "hhn")
fascr.hr <- fit.ascr(fcapt.nullSS, ftraps, fmask, detfn = "hr")

# frog models without auxiliary toa data ----------------------------------
fcapt.nullTOA <- list(bincapt = fcapt$bincapt, ss = fcapt$ss)
fascr.noTOA <- fit.ascr(fcapt.nullTOA, ftraps, fmask, detfn = "ss", ss.opts = list(cutoff = 140, ss.link = "identity"))

fcapt.nullAUX <- list(bincapt = fcapt$bincapt)
fascr.hn.noTOA <- fit.ascr(fcapt.nullAUX, ftraps, fmask, detfn = "hn")
fascr.hhn.noTOA <- fit.ascr(fcapt.nullAUX, ftraps, fmask, detfn = "hhn")
fascr.hr.noTOA <- fit.ascr(fcapt.nullAUX, ftraps, fmask, detfn = "hr", sv = list(D = 1000)) ## Note that we fix it

# Save the models ---------------------------------------------------------
sim.theta <- list(
  fascr = list(theta = get.par(fascr, pars = "fitted", as.list = TRUE), detfn = "ss", ss.opts = fascr$args$ss.opts, infotypes = "toa")
)
saveRDS(sim.theta, file = "./code/miss-v2/model_fit.RDS")

# -------------------------------------------------------------------------

fascr.Dc = data.frame(point = c(coef(fascr.hn)["D"], coef(fascr.hr)["D"], coef(fascr.hhn)["D"], coef(fascr)["D"]),
                      lower = c(confint(fascr.hn)["D", 1], confint(fascr.hr)["D", 1], confint(fascr.hhn)["D", 1], confint(fascr)["D", 1]),
                      upper = c(confint(fascr.hn)["D", 2], confint(fascr.hr)["D", 2], confint(fascr.hhn)["D", 2], confint(fascr)["D", 2])) / 25

fascr.Dc.2 = data.frame(point = c(coef(fascr.hn.noTOA)["D"], coef(fascr.hr.noTOA)["D"], coef(fascr.hhn.noTOA)["D"], coef(fascr.noTOA)["D"]),
                        lower = c(confint(fascr.hn.noTOA)["D", 1], confint(fascr.hr.noTOA)["D", 1], confint(fascr.hhn.noTOA)["D", 1], confint(fascr.noTOA)["D", 1]),
                        upper = c(confint(fascr.hn.noTOA)["D", 2], confint(fascr.hr.noTOA)["D", 2], confint(fascr.hhn.noTOA)["D", 2], confint(fascr.noTOA)["D", 2])) / 25

# -------------------------------------------------------------------------

pdf("./figures/05-figureN1.pdf", height = 5, width = 8)
palette("Dark 2")
par(mai = c(1, 0.55, 0.15, 0.15), xaxs = "i", oma = rep(0.3, 4))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(0.5, 4.5), ylim = range(fascr.Dc, fascr.Dc.2))

abline(h = seq(0, 200, 10), col = "grey90", lty = 1)

points(x = 1:4 - 0.25, y = fascr.Dc$point, cex = 2, pch = 18, col = c(4,6,3,1))
segments(x0 = 1:4 - 0.25, y0 = fascr.Dc[, 2], y1 = fascr.Dc[, 3], lwd = 2, col = c(4,6,3,1))
segments(x0 = 0.9 + 0:3 - 0.25, y0 = fascr.Dc[, 2], x1 = 1.1 + 0:3 - 0.25, lwd = 2, col = c(4,6,3,1))
segments(x0 = 0.9 + 0:3 - 0.25, y0 = fascr.Dc[, 3], x1 = 1.1 + 0:3 - 0.25, lwd = 2, col = c(4,6,3,1))
segments(x0 = 1:4 - 0.25, y0 = fascr.Dc[, 2], y1 = fascr.Dc[, 3], lwd = 2, col = c(4,6,3,1))

points(x = 1:4 + 0.25, y = fascr.Dc.2$point, cex = 2, pch = 18, col = c(4,6,3,1))
segments(x0 = 1:4 + 0.25, y0 = fascr.Dc.2[, 2], y1 = fascr.Dc.2[, 3], lwd = 2, col = c(4,6,3,1), lty = 6)
segments(x0 = 0.9 + 0:3 + 0.25, y0 = fascr.Dc.2[, 2], x1 = 1.1 + 0:3 + 0.25, lwd = 2, col = c(4,6,3,1))
segments(x0 = 0.9 + 0:3 + 0.25, y0 = fascr.Dc.2[, 3], x1 = 1.1 + 0:3 + 0.25, lwd = 2, col = c(4,6,3,1))
segments(x0 = 1:4 + 0.25, y0 = fascr.Dc.2[, 2], y1 = fascr.Dc.2[, 3], lwd = 2, col = c(4,6,3,1), lty = 6)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0), at = c(1:4), labels = c("Halfnormal", "Hazard rate", "Hazard-halfnormal", "Signal strength")); mtext("Detection function", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0), at = seq(0, 100, 10));

mtext(expression(bold("Call density p/ha p/s")), side = 2, line = 1.6, font = 2, cex = 1.1)

box()

par(fig = c(0, 1, 0, 0.975), oma = c(0, 0, 0, 0), mar = c(0.25, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', legend = c("with Times of Arrival", "without Times of Arrival"), lty = c(1, 6), col = "black", lwd = 2,
       horiz = TRUE, cex = 0.9)

dev.off()

# -------------------------------------------------------------------------

pdf("./figures/05-figureN2.pdf", height = 6.5, width = 12)
palette("Dark 2")

layout(matrix(c(1, 2), 1, 2))

par(mai = c(1, 0.55, 0.25, 0.15), xaxs = "i", oma = rep(0.3, 4))

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(0, 32.5), ylim = c(0, 1), main = "with Times of Arrival")

abline(h = seq(0, 1, 0.2), v = seq(0, 30, 5), col = "grey90", lty = 1)

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

show.detfn(fascr.hn, add = TRUE, col = 4, lwd = 2)
show.detfn(fascr.hr, add = TRUE, col = 6, lwd = 2)
show.detfn(fascr.hhn, add = TRUE, col = 3, lwd = 2)
show.detfn(fascr, add = TRUE, col = 1, lwd = 2)

box()

# legend("topright", bg = "white", legend = c("Halfnormal", "Hazard rate", "Hazard-halfnormal", "Signal strength"), lty = 1, col = 1:4, lwd = 2)

#

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = c(0, 32.5), ylim = c(0, 1), main = "without Times of Arrival")

abline(h = seq(0, 1, 0.2), v = seq(0, 30, 5), col = "grey90", lty = 1)


axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.7, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.3, 0)); mtext("Probability of detection", side = 2, line = 1.6, font = 2, cex = 1.1)

show.detfn(fascr.hn.noTOA, add = TRUE, col = 4, lwd = 2, lty = 6)
show.detfn(fascr.hr.noTOA, add = TRUE, col = 6, lwd = 2, lty = 6)
show.detfn(fascr.hhn.noTOA, add = TRUE, col = 3, lwd = 2, lty = 6)
show.detfn(fascr.noTOA, add = TRUE, col = 1, lwd = 2, lty = 6)

box()

# legend("topright", bg = "white", legend = c("Halfnormal", "Hazard rate", "Hazard-halfnormal", "Signal strength"), lty = 6, col = 1:4, lwd = 2)

par(fig = c(0, 1, 0, 0.975), oma = c(0, 0, 0, 0), mar = c(0.25, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom', bg = "white", legend = c("Halfnormal", "Hazard rate", "Hazard-halfnormal", "Signal strength"), lty = 1, col = c(4,6,3,1), lwd = 2,
       horiz = TRUE)

dev.off()

c(hn = AIC(fascr.hn.noTOA), hr = AIC(fascr.hr.noTOA), hhn = AIC(fascr.hhn.noTOA)) |> sort()

c(hn = AIC(fascr.hn), hr =  AIC(fascr.hr), hhn = AIC(fascr.hhn)) |> sort()
