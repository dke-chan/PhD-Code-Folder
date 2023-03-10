source("./code/gibbons/helper-Scripts/simul.r")

library(dplyr)

# Data
set.seed(2739282)
posts.input <- expand.grid(x = c(-1000, 0, 1000), y = c(0))
# The.detfn.theta <- c(lambda0 = 5.75, sigma = 500)
The.detfn.theta <- c(lambda0 = 3, sigma = 550)
The.move.theta <- diag(550, 2) ^ 2 |> as.vector()
survey.length <- 5
simul.frame <- simul_data(survey.region = data.frame(x = 16000 * c(-1, 1), y = 15000 * c(-1, 1)),
                          listening.post = posts.input, # movement.theta = rep(0, 4),
                          detfn.theta = The.detfn.theta,
                          movement.theta = The.move.theta,
                          # call.Process = "Beronulli", call.Process.Theta = c(survey.length, 2/3))
                          call.Process = "Binomial", call.Process.Theta = c(survey.length, 2/3))
                          # call.Process = "Poisson", call.Process.Theta = c(1.25))

# Construct the required matrices -----------------------------------------

posts <- as.matrix(posts.input)
mask <- ascr::create.mask(posts, buffer = 10000)
mask.dists <- ascr::distances(mask, posts)
expected.bearing <- t(ascr::bearings(posts, mask))

# Create the capture history and their bearings from `simul.frame`
bin.capt <- simul.frame %>%
  mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
  filter(detected > 0) %>%
  select(starts_with("detect.")) %>%
  as.matrix()

bearings <- simul.frame %>%
  mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
  filter(detected > 0) %>%
  select(starts_with("ob.bearing")) %>%
  as.matrix()

# Pull the group identifiers out
animal.id <- simul.frame %>%
  mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>%
  filter(detected > 0) %>%
  pull(groupID) %>%
  factor() %>%
  as.numeric() %>%
  {. - 1}

# Convert the observed bearings to be 0s for non-detections
bearings[, 1:nrow(posts)] <- bearings[, 1:nrow(posts)] * bin.capt[, 1:nrow(posts)]

# Build local versions of the Rcpp functions ------------------------------

library(Rcpp); library(RcppArmadillo)
sourceCpp("./code/gibbons/cpp-Files/ascr.cpp")
sourceCpp("./code/gibbons/cpp-Files/ascrStack.cpp"); sourceCpp("./code/gibbons/cpp-Files/ascrDisperse.cpp")

# Standard ASCR with no animal identities ---------------------------------

std.ascr.start <- c(log(0.75 * (2/3 * survey.length) / 100), log(3), log(1250), log(36.22))
std.ascr <- optim(std.ascr.start, fit_ascr, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask, traps = posts,
                  observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
                  mask_dists = mask.dists, pixel_area = attr(mask, "area"), hessian = TRUE, method = "L-BFGS-B")

std.ascr.res <- cbind(theta.hat = exp(std.ascr$par),
                      true.theta = c(D_call = 0.75 * (2/3 * survey.length), lambda0 = 3, sigma = 1250, kappa = 36.22),
                      wald.lwr = exp(std.ascr$par - qnorm(0.975) * sqrt(diag(solve(std.ascr$hessian)))),
                      wald.upr = exp(std.ascr$par + qnorm(0.975) * sqrt(diag(solve(std.ascr$hessian)))))
std.ascr.res[1, -2] <- std.ascr.res[1, -2] * 100

# ASCR with animal identities & stacking ~ Stevenson et al ----------------

stack.ascr.start <- c(log(0.75 / 100), log(3), log(1250), qlogis(2/3), log(36.22))
stack.ascr <- optim(stack.ascr.start, fit_ascrStack, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask,
                    traps = posts, observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
                    mask_dists = mask.dists, pixel_area = attr(mask, "area"), ID = animal.id,
                    survey_length = survey.length, bernoulli = TRUE, hessian = TRUE, method = "L-BFGS-B")

stack.ascr.hess <- try(diag(solve(stack.ascr$hessian)))
if (class(stack.ascr.hess) == "try-error") {
  stack.ascr.hess <- c(diag(solve(stack.ascr$hessian[1:4, 1:4])), 0)
}
stack.ascr.res <- cbind(theta.hat = stack.ascr$par,
                        true.theta = c(D = 0.75, lambda0 = 3, sigma = 1250, pBernoulli = 2/3, kappa = 36.22),
                        wald.lwr = stack.ascr$par - qnorm(0.975) * sqrt(stack.ascr.hess),
                        wald.upr = stack.ascr$par + qnorm(0.975) * sqrt(stack.ascr.hess))
stack.ascr.res[-4, -2] <- exp(stack.ascr.res[-4, -2])
stack.ascr.res[4, -2] <- plogis(stack.ascr.res[4, -2])
stack.ascr.res[1, -2] <- stack.ascr.res[1, -2] * 100

# ASCR with animal identities and dispersal -------------------------------
disperse.ascr.start <- c(log(0.75 / 100), log(3), log(1250), qlogis(2/3), log(1000), log(36.22))
disperse.ascr <- optim(disperse.ascr.start, fit_ascrDisperse, bin_capt = as.matrix(bin.capt[, 1:nrow(posts)]), mask = mask,
                       traps = posts, observed_bearings = as.matrix(bearings[, 1:nrow(posts)]), expected_bearing = expected.bearing,
                       mask_dists = mask.dists, pixel_area = attr(mask, "area"), ID = animal.id, survey_length = survey.length,
                       hessian = TRUE, method = "L-BFGS-B", trace = TRUE)

disperse.ascr.res <- cbind(theta.hat = disperse.ascr$par,
                           true.theta = c(D = 0.75, lambda0 = 3, sigma = 1250, pBernoulli = 2/3, sigmaMove = 1000, kappa = 36.22),
                           wald.lwr = disperse.ascr$par - qnorm(0.975) * sqrt(diag(solve(disperse.ascr$hessian))),
                           wald.upr = disperse.ascr$par + qnorm(0.975) * sqrt(diag(solve(disperse.ascr$hessian))))
disperse.ascr.res[-4, -2] <- exp(disperse.ascr.res[-4, -2])
disperse.ascr.res[4, -2] <- plogis(disperse.ascr.res[4, -2])
disperse.ascr.res[1, -2] <- disperse.ascr.res[1, -2] * 100

# nll: nan D: 0.00508867 lambda0: 2.58339 sigma: 479.631 pBinomial: 0.872152 sigmaMove: 41.8682 kappa: 8.52103

# -------------------------------------------------------------------------

vcov.hat <- solve(disperse.ascr$hessian)
G = diag(length(disperse.ascr$par)) * c(exp(disperse.ascr$par[1:3]),  # deriv of exp()
                                        dlogis(disperse.ascr$par[4]), # deriv of plogis()
                                        exp(disperse.ascr$par[5:6]))  # deriv of exp()
se = sqrt(diag(G %*% vcov.hat %*% t(G))) ## For the back-transformed

round(std.ascr.res, 2)
round(stack.ascr.res, 2)
round(disperse.ascr.res, 2)

matrix(survey.length * c(plogis(disperse.ascr$par[4]), exp(disperse.ascr$par[1])), 1, 2) %*% (G %*% vcov.hat %*% t(G))[c(1, 4), c(1, 4)] %*%  t(matrix(survey.length * c(plogis(disperse.ascr$par[4]), exp(disperse.ascr$par[1])), 1, 2)) |>
  sqrt()

# -------------------------------------------------------------------------

G = c(survey.length * exp(disperse.ascr$par[1]) * plogis(disperse.ascr$par[4]),
      survey.length * exp(disperse.ascr$par[1]) * exp(disperse.ascr$par[4]) / (exp(disperse.ascr$par[4]) + 1)^2) |>
  matrix(1, 2)
sqrt(G %*% vcov.hat[c(1, 4), c(1, 4)] %*% t(G)) |>
  as.numeric() ## For call density



# msm::deltamethod(list(~ exp(x1), ~ exp(x4) / (exp(x4) + 1)), disperse.ascr$par, vcov.hat, ses = FALSE)

# G[c(1, 4), c(1, 4)] %*% vcov.hat[c(1, 4), c(1, 4)] %*% t(G[c(1, 4), c(1, 4)])
exp(disperse.ascr$par[1]) * dlogis(disperse.ascr$par[4]) * vcov.hat[1, 4]


toReport <- rbind(
  ascrDisperse = c(disperse.ascr.res[1] * disperse.ascr.res[4] * survey.length, disperse.ascr.res[1] * disperse.ascr.res[4] * survey.length + se * c(-1, 1) * qnorm(0.975)),
  ascr = c(std.ascr.res[1, 1], std.ascr.res[1, 3:4]))

toReport


