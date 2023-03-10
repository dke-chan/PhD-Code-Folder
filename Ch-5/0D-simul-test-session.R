source("./code/gibbons/helper-Scripts/simul.r")

library(dplyr); library(ascr)

# Data
set.seed(2739282)
posts.input <- list(expand.grid(x = c(-11000, -10000, -9000), y = 7500),
                    expand.grid(x = c(11000, 10000, 9000), y = 7500),
                    expand.grid(x = c(-1000, 0, 1000), y = c(0)),
                    expand.grid(x = c(-11000, -10000, -9000), y = -7500),
                    expand.grid(x = c(11000, 10000, 9000), y = -7500))
The.detfn.theta <- c(lambda0 = 3, sigma = 550)
The.move.theta <- diag(1000, 2) ^ 2 |> as.vector()
survey.length <- 3
simul.frame <- simul_session(n.session = 5, listening.post.session = posts.input,
                             survey.region = data.frame(x = 16000 * c(-1, 1), y = 15000 * c(-1, 1)),
                             # movement.theta = rep(0, 4),
                             detfn.theta = The.detfn.theta,
                             movement.theta = The.move.theta,
                             call.Process = "Binomial", call.Process.Theta = c(survey.length, 2/3))

# -------------------------------------------------------------------------

posts.session <- lapply(posts.input, as.matrix)
mask.session <- lapply(posts.session, create.mask, buffer = 7500)
mask.dist.session <- lapply(1:length(mask.session), \(x) {ascr::distances(mask.session[[x]], posts.session[[x]])})
expected.bearing.session <- lapply(1:length(mask.session), \(x) {t(ascr::bearings(posts.session[[x]], mask.session[[x]]))})
pixel.size.session <- sapply(mask.session, attr, "area")

bincapt.session <- lapply(simul.frame, \(x) {
  y <- x[, grep("detect\\.", colnames(x))]
  y[which(rowSums(y) > 0), ] |> as.matrix()
})

bearing.session <- lapply(simul.frame, \(x) {
  y <- x[, grep("ob.bearing\\.", colnames(x))]
  z <- x[, grep("detect\\.", colnames(x))]
  as.matrix(y[which(rowSums(z) > 0), ] * z[which(rowSums(z) > 0), ])
})

animal.id.session <- lapply(simul.frame, \(x) {
  y <- x[, grep("detect\\.", colnames(x))]
  z <- x$groupID[which(rowSums(y) > 0)] |> factor() |> as.numeric()
  z - 1
})

# -------------------------------------------------------------------------
sourceCpp("./code/gibbons/cpp-Files/ascrStack.cpp"); sourceCpp("./code/gibbons/cpp-Files/ascrDisperse.cpp")

# -------------------------------------------------------------------------

capt.calls <- lapply(1:length(bincapt.session), \(x) {
  list(bincapt = bincapt.session[[x]], bearing = bearing.session[[x]])
})

benchmark.fit <- fit.ascr(capt.calls, traps = posts.session, detfn = "hhn", mask = mask.session)

# -------------------------------------------------------------------------

disperse.start <- c(log(0.75 / 100), log(3), log(550), qlogis(2/3), log(1000), log(36.22))
disperse.model <- optim(disperse.start, fit_ascrDisperseSession, bin_capt = bincapt.session, mask = mask.session,
                        traps = posts.session, mask_dists = mask.dist.session, observed_bearings = bearing.session,
                        expected_bearings = expected.bearing.session, pixel_area = pixel.size.session,
                        ID = animal.id.session, survey_length = rep(survey.length, length(bincapt.session)), threshold = 1e-4,
                        trace = TRUE, hessian = TRUE, method = "Nelder-Mead")

vcov.hat <- solve(disperse.model$hessian)
G = diag(length(disperse.model$par)) * c(exp(disperse.model$par[1:3]),
                                         plogis(disperse.model$par[4]),
                                         exp(disperse.model$par[5:6]))
se = sqrt(diag(G %*% vcov.hat %*% t(G))) ## For the back-transformed


G = c(survey.length * exp(disperse.model$par[1]) * plogis(disperse.model$par[4]),
      rep(0, 2),
      survey.length * exp(disperse.model$par[1]) * exp(disperse.model$par[4]) / (exp(disperse.model$par[4]) + 1)^2,
      rep(0, 2))
se = sqrt(t(G) %*% vcov.hat %*% G) |>
  as.numeric() ## For call density

disperse.ascr.res <- cbind(theta.hat = disperse.model$par,
                           true.theta = c(D = 0.75, lambda0 = 3, sigma = 550, pBernoulli = 2/3, sigmaMove = 1000, kappa = 36.22),
                           wald.lwr = disperse.model$par - qnorm(0.975) * sqrt(diag(solve(disperse.model$hessian))),
                           wald.upr = disperse.model$par + qnorm(0.975) * sqrt(diag(solve(disperse.model$hessian))))
disperse.ascr.res[-4, -2] <- exp(disperse.ascr.res[-4, -2])
disperse.ascr.res[4, -2] <- plogis(disperse.ascr.res[4, -2])
disperse.ascr.res[1, -2] <- disperse.ascr.res[1, -2] * 100


# -------------------------------------------------------------------------

cbind(Est. = 100 * coef(benchmark.fit), Truth = c(D = 0.75 * 2/3 * survey.length, lambda0 = 3, sigma = 550, kappa = 36.22), 100 * confint(benchmark.fit))
disperse.ascr.res

G = c(survey.length * exp(disperse.model$par[1]) * plogis(disperse.model$par[4]),
      rep(0, 2),
      survey.length * exp(disperse.model$par[1]) * exp(disperse.model$par[4]) / (exp(disperse.model$par[4]) + 1)^2,
      rep(0, 2))
se = sqrt(t(G) %*% vcov.hat %*% G) |>
  as.numeric() ## For call density

toReport <- rbind(
  ascrDisperse = c(disperse.ascr.res[1] * disperse.ascr.res[4] * survey.length, disperse.ascr.res[1] * disperse.ascr.res[4] * survey.length + se * c(-1, 1) * qnorm(0.975)),
  ascr = 100 * c(coef(benchmark.fit)[1], confint(benchmark.fit)[1, ]))

toReport

plot(0, xlim = c(0, 3), ylim = c(0, 6), type = "n")
points(x = 1:2, y = toReport[, 1], pch = 16)
points(x = rep(1:2, 2), y = c(toReport[, 2], toReport[, 3]), pch = 16)
