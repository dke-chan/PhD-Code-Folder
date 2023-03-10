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
    g0 <- exp(theta[1])/(exp(theta[1]) + 1); sigma <- exp(theta[2]); zeta <- exp(theta[3])
  } else {
    g0 <- theta[1]; sigma <- theta[2]; zeta <- theta[3]
  }

  return (g0 * (1 - exp(-(d / sigma)^(-zeta))))
}

threshold <- function (d, theta, link = TRUE) {
  if (link) {
    shape <- theta[1]; scale <- exp(theta[2])
  } else {
    shape <- theta[1]; scale <- theta[2]
  }

  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

  return (0.5 - 0.5 * erf(d/scale - shape))
}

# -------------------------------------------------------------------------

bernoulliLL <- function (theta, y, x, probFunc) {
  psi <- log(probFunc(x, theta) + .Machine$double.eps) - log(1 - probFunc(x, theta) + .Machine$double.eps);
  return (-sum((y * psi - log(1 + exp(psi) + .Machine$double.eps))))
}

# -------------------------------------------------------------------------

modelResPrinter <- function (x, w) {
  if (inherits(x, "character")) {
    y <- get(x)
  } else {
    y <- x
  }

  a <- confint(y)
  colnames(a) <- c("lwr", "upr")

  b <- cbind(est = coef(y), a) |>
    as.data.frame()

  b$detfn <- y$args$detfn

  b$est.name <- names(coef(y))

  b$spec <- w

  b$se <- y$se[b$est.name]

  rownames(b) <- NULL

  return (b)
}


# -------------------------------------------------------------------------

simNewCaptures <- function(detfn, theta, kappa) {
  detfn.f <- switch(detfn, hn = halfnormal, hr = hazardRate,
                    th = threshold, hhn = hazardHalfnormal)

  capt.starAll <- sapply(full.distance.mat, \(x) {
    rbinom(1, 1, detfn.f(x * 1000, theta, link = FALSE))
  }) |> matrix(nrow = 30, ncol = 4)

  bearing.starAll <- sapply(full.true.bearing.mat, \(x) {
    rvm(1, x, k = kappa)
  }) |> matrix(nrow = 30, ncol = 4)

  whichRows <- which(rowSums(capt.starAll) > 0, arr.ind = TRUE)

  subset(seaTrial.df, PlaybackNumber == 1)[, c("buoy.rel.x", "buoy.rel.y")] |>
    as.matrix()

  return(list(
    capt.obj = lapply(whichRows, \(x) {
      list(
        bincapt = capt.starAll[x, , drop = FALSE],
        bearing = bearing.starAll[x, , drop = FALSE],
        mrds = full.boat.ave.rel.mat[x, , drop = FALSE]
      )
    }),
    traps.obj = lapply(whichRows, \(x) {
      subset(seaTrial.df, PlaybackNumber == x)[, c("buoy.rel.x", "buoy.rel.y")] |>
        as.matrix()
    }),
    mask.obj = lapply(vector("list", length(whichRows)), function(x) mask.obj),
    capt.obj.full = capt.starAll))
}

# -------------------------------------------------------------------------

euclideanSq <- function (d, theta.1, theta.2, detfn, link = FALSE) {
  (detfn(d, theta.1, link = FALSE) - detfn(d, theta.2, link = FALSE))^2
}

euclideanNorm <- function(theta.1, theta.2, detfn, link = FALSE) {
  return (pracma::quadinf(euclideanSq, 0, Inf, theta.1 = theta.1, theta.2 = theta.2, detfn = detfn, link = link)$Q^(1/2))
}

