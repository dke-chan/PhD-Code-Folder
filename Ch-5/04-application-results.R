library(readxl); library(ascr); library(secr); library(Rcpp); library(RcppArmadillo)

load("code/gibbons/01-data-prep.RData")
load("code/gibbons/03-fit-models.RData")

sourceCpp("code/gibbons/cpp-Files/ascrStack.cpp"); sourceCpp("code/gibbons/cpp-Files/ascrDisperse.cpp")

# ascr --------------------------------------------------------------------
ascr.results <- cbind(
  Theta = summary(benchmark.fit)$coefs,
  SE = summary(benchmark.fit)$coefs.se,
  CV = 100 * summary(benchmark.fit)$coefs.se / summary(benchmark.fit)$coefs,
  lwr.Wald = confint(benchmark.fit, linked = TRUE)[, 1],
  upr.Wald = confint(benchmark.fit, linked = TRUE)[, 2]
  ) |>
  rbind(c(mean(benchmark.derived[grep("esa", names(benchmark.derived))]) / 100, NA, NA, NA, NA))
rownames(ascr.results)[5] <- "esa_bar"
ascr.results[1, c(1, 2, 4, 5)] <- 100 * ascr.results[1, c(1, 2, 4, 5)] / 3

# ascrStack ---------------------------------------------------------------
## estimated
stack.theta <- c(100 * exp(stack.model$par[1]), exp(stack.model$par[2:3]), plogis(stack.model$par[4]), exp(stack.model$par[5]))

G <- diag(c(100 * exp(stack.model$par[1]), exp(stack.model$par[2:3]), dlogis(stack.model$par[4]), exp(stack.model$par[5])))

est.vcov.hat <- G %*% solve(stack.model$hessian) %*% t(G)

stack.se <- sqrt(diag(est.vcov.hat))

walds <- stack.model$par + cbind(rep(-1, length(stack.theta)), rep(1, length(stack.theta))) * qnorm(0.975) * sqrt(diag(solve(stack.model$hessian)))
walds[-4, ] <- exp(walds[-4, ])
walds[4, ] <- plogis(walds[4, ])
walds[1, ] <- 100 * walds[1, ]

## derived
Dc.hat <- survey.length * stack.theta[1] * stack.theta[4] / survey.length
G.Dc <- matrix(100 * survey.length / survey.length *
                 c(exp(stack.model$par[1]) * plogis(stack.model$par[4]), 0, 0,
                   exp(stack.model$par[1]) * exp(stack.model$par[4]) / (exp(stack.model$par[4]) + 1)^2, 0),
               1, 5)
Dc.se <- sqrt(G.Dc %*% solve(stack.model$hessian) %*% t(G.Dc))
Dc.walds <- Dc.hat + c(-1, 1) * qnorm(0.975) * as.numeric(Dc.se)

ascrStack.esa_bar <- sapply(1:length(posts.full), \(x) {
  dists <- distances(survey.session[[x]], posts.full[[x]])

  call.probs <- 1 - exp(-exp(stack.model$par[2])*exp(-dists^2/(2*exp(stack.model$par[3])^2)))

  p_c <- 1 - apply(1 - call.probs, 1, prod)

  ind.probs <- 1 - (1 - plogis(stack.model$par[4]) * p_c)^survey.length

  sum(ind.probs * attr(survey.session[[x]], "area"))
}) |> mean() / 100

ascrStack.results <- cbind(Theta = stack.theta, SE = stack.se, CV = 100 * stack.se / stack.theta, lwr.Wald = walds[, 1], upr.Wald = walds[, 2]) |>
  rbind(cbind(Theta = c(Dc.hat), SE = c(Dc.se), CV = 100 * c(Dc.se) / c(Dc.hat), lwr.Wald = c(Dc.walds[1]), upr.Wald = c(Dc.walds[2]))) |>
  rbind(c(ascrStack.esa_bar, NA, NA, NA, NA))
rownames(ascrStack.results) <- c("D", "lambda0", "sigma", "pb", "kappa", "Dc", "esa_bar")

# ascrDisperse ------------------------------------------------------------
## estimated
disperse.theta <- c(100 * exp(disperse.model$par[1]), exp(disperse.model$par[2:3]), plogis(disperse.model$par[4]), exp(disperse.model$par[5:6]))

G <- diag(c(100 * exp(disperse.model$par[1]), exp(disperse.model$par[2:3]), dlogis(disperse.model$par[4]), exp(disperse.model$par[5:6])))
est.vcov.hat <- G %*% solve(disperse.model$hessian) %*% t(G)

disperse.se <- sqrt(diag(est.vcov.hat))

walds <- disperse.model$par + cbind(rep(-1, length(disperse.theta)), rep(1, length(disperse.theta))) * qnorm(0.975) * sqrt(diag(solve(disperse.model$hessian)))
walds[-4, ] <- exp(walds[-4, ])
walds[4, ] <- plogis(walds[4, ])
walds[1, ] <- 100 * walds[1, ]

## derived
Dc.hat <- survey.length / survey.length * disperse.theta[1] * disperse.theta[4]
G.Dc <- matrix(100 * survey.length  / survey.length *
                 c(exp(disperse.model$par[1]) * plogis(disperse.model$par[4]), 0, 0,
                   exp(disperse.model$par[1]) * exp(disperse.model$par[4]) / (exp(disperse.model$par[4]) + 1)^2, 0, 0),
               1, 6)
Dc.se <- sqrt(G.Dc %*% solve(disperse.model$hessian) %*% t(G.Dc))

Dc.walds <- Dc.hat + c(-1, 1) * qnorm(0.975) * as.numeric(Dc.se)

ra.hat <- pi * (exp(disperse.model$par[5]) * 2.4477)^2 / 10000
G.ra <- matrix(pi / 10000 * c(0, 0, 0, 0, 2.4477 * 2.4477 * 2 * exp(disperse.model$par[5])^2, 0), 1, 6)
ra.se <- sqrt(G.ra %*% solve(disperse.model$hessian)%*% t(G.ra))

ra.walds <- disperse.model$par[5] + c(-1, 1) * qnorm(0.975) * sqrt(diag(solve(disperse.model$hessian)))[5]
ra.walds <- (exp(ra.walds) * 2.4477)^2 * pi / 10000

## call dispersal average esa_c in km^2
ascr.Disperse.esac_bar <- sapply(1:length(posts.full), \(x) {
  dists <- distances(survey.session[[x]], posts.full[[x]])

  call.probs <- 1 - exp(-exp(disperse.model$par[2])*exp(-dists^2/(2*exp(disperse.model$par[3])^2)))

  p_c <- 1 - apply(1 - call.probs, 1, prod)

  sum(p_c * attr(survey.session[[x]], "area"))
}) |> mean() / 100

## call dispersal average esa_a in km^2
ascr.Disperse.esa_bar <- sapply(1:length(posts.full), \(x) {
  dists <- distances(survey.session[[x]], posts.full[[x]])

  call.probs <- 1 - exp(-exp(disperse.model$par[2])*exp(-dists^2/(2*exp(disperse.model$par[3])^2)))

  p_c <- 1 - apply(1 - call.probs, 1, prod)

  biNormConstraint <- qnorm((2 - 1e-4) / 2) * exp(disperse.model$par[5])

  m2_pixel_area <- attr(survey.session[[x]], "area") * 10000
  m2_pixel_edge <- sqrt(m2_pixel_area)

  p_actCen <- apply(survey.session[[x]], 1, \(a) {
    inner.dists <- distances(survey.session[[x]], matrix(a, 1, 2))

    indices <- which(inner.dists < biNormConstraint)

    if (biNormConstraint > m2_pixel_edge) {
      fs <- m2_pixel_area * p_c[indices] * mvtnorm::dmvnorm(survey.session[[x]][indices, ], as.numeric(a), diag(exp(disperse.model$par[5]), 2)^2)
    } else {
      fs <- p_c[indices]
    }

    return (sum(fs))
  })

  ind.probs <- 1 - (1 - plogis(disperse.model$par[4]) * p_actCen)^survey.length

  sum(ind.probs * attr(survey.session[[x]], "area"))
}) |> mean() / 100

ascrDisperse.results <- cbind(Theta = disperse.theta, SE = disperse.se, CV = 100 * disperse.se / disperse.theta, lwr.Wald = walds[, 1], upr.Wald = walds[, 2]) |>
  rbind(cbind(Theta = c(Dc.hat, ra.hat), SE = c(Dc.se, ra.se), CV = 100 * c(Dc.se, ra.se) / c(Dc.hat, ra.hat), wr.Wald = c(Dc.walds[1], ra.walds[1]), upr.Wald = c(Dc.walds[2], ra.walds[2]))) |>
  rbind(cbind(Theta = c(ascr.Disperse.esac_bar, ascr.Disperse.esa_bar), NA, NA, NA, NA))
rownames(ascrDisperse.results) <- c("D", "lambda0", "sigma", "pb", "sigmaMove", "kappa", "Dc", "ra", "ebarc", "ebara")

# Figure 4.2 --------------------------------------------------------------

dists <- 0:5250

benchmark.gc <- ascr:::calc.hhn(dists, list(lambda0 = coef(benchmark.fit)[2], sigma = coef(benchmark.fit)[3]))
stack.gc <- ascr:::calc.hhn(dists, list(lambda0 = stack.theta[2], sigma = stack.theta[3]))
disperse.gc <- ascr:::calc.hhn(dists, list(lambda0 = disperse.theta[2], sigma = disperse.theta[3]))

# stack.ga <- 1 - (1 - stack.theta[4] * stack.gc)^survey.length

model.cols <- RColorBrewer::brewer.pal(3, "Dark2")

pdf("./figures/04-figure02.pdf", height = 6, width = 6)
par(mai = c(0.45, 0.5, 0.05, 0.05), xaxs = "i")
plot(0, type = "l", axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1), xlim = c(0, 3250))
abline(h = c(0, 1), col = "grey75")

lines(x = dists, y = benchmark.gc, col = model.cols[1], lwd = 2)
lines(x = dists, y = stack.gc, col = model.cols[2], lwd = 2)
lines(x = dists, y = disperse.gc, col = model.cols[3], lwd = 2)

# lines(x = dists, y = stack.ga, col = model.cols[2], lwd = 2, lty = 2)

legend("topright", bg = "white", lty = 1, legend = c("Standard", "Call Stacking", "Call Dispersal"), lwd = 2, col = model.cols[1:3])

axis(1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Distance (m)", side = 1, line = 1.25, font = 2, cex = 0.75)
axis(2, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.4, 0)); mtext("Probability of detection", side = 2, line = 1.55, font = 2, cex = 0.8)

box()
dev.off()

# -------------------------------------------------------------------------

ascr.results
ascrStack.results
ascrDisperse.results
