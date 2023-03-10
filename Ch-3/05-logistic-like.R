load("./data/sonobuoys.RData")
source("./code/sonobuoys/00-aux-funcs.R")

# -------------------------------------------------------------------------
bernoulli.trial <- as.vector(full.capt.hist.mat)
bernoulli.distance <- as.vector(seaTrial.df$Distance * 1000)

##
logi.hn <- optim(c(1, log(1000)), bernoulliLL, probFunc = halfnormal, y = bernoulli.trial, x = bernoulli.distance, hessian = TRUE, method = "L-BFGS-B")
logi.hr <- optim(c(1, log(1000), log(10)), bernoulliLL, probFunc = hazardRate, y = bernoulli.trial, x = bernoulli.distance, hessian = TRUE)
logi.th <- optim(c(1, log(1000)), bernoulliLL, probFunc = threshold, y = bernoulli.trial, x = bernoulli.distance, hessian = TRUE, method = "L-BFGS-B")
logi.hhn <- optim(c(1, log(1000)), bernoulliLL, probFunc = hazardHalfnormal, y = bernoulli.trial, x = bernoulli.distance, hessian = TRUE, method = "L-BFGS-B")

# Compiling for 07-inference.r --------------------------------------------

link.ests <- lapply(ls(pattern = c("logi\\.")), \(x) {
  y <- get(x)
  hess <- solve(y$hessian)
  whichFun <- gsub("logi\\.(*)", "\\1", x)

  z <- data.frame(est = y$par,
             lwr = 0,
             upr = 0,
             detfn = whichFun,
             est.name = switch(whichFun, hn = c("g0", "sigma"), hhn = c("lambda0", "sigma"), hr = c("g0", "sigma", "z"), th = c("shape", "scale") ),
             spec = "drift.unbatched")

  to.exp <- which(z$est.name %in% c("sigma", "lambda0", "scale", "z"))
  to.plogis <- which(z$est.name %in% c("g0"))

  z$est[to.exp] <- exp(y$par[to.exp])
  z$est[to.plogis] <- plogis(y$par[to.plogis])

  G.mult <- y$par
  G.mult[to.exp] <- exp(G.mult[to.exp])
  G.mult[to.plogis] <- dlogis(G.mult[to.plogis])

  G <- diag(length(y$par)) * G.mult

  se <- sqrt(diag(G %*% hess %*% t(G)))

  z$lwr <- z$est - se * qnorm(0.975)
  z$upr <- z$est + se * qnorm(0.975)

  z$se <- se

  return (z)
}) |> do.call(what = "rbind")

saveRDS(link.ests, "./code/sonobuoys/logistic-like-models.RDS")
