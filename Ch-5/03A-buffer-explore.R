library(readxl); library(ascr); library(secr); library(Rcpp); library(RcppArmadillo)
load("code/gibbons/01-data-prep.RData")

# -------------------------------------------------------------------------

disperse.start <- c(log(0.005), log(2), log(650), qlogis(0.5), log(550), log(8.28))

buffer.choice <- seq(1000, 12500, by = 500)
# spacing.choice <- seq(150, 300, by = 50)
maskSpecs <- expand.grid(buffer = buffer.choice) #, spacing = spacing.choice)
maskSpecs.List <- split(maskSpecs, seq(nrow(maskSpecs)))

# -------------------------------------------------------------------------

if (!file.exists("code/gibbons/03A-buffer-explore.rds")) {
  library(parallel)
  cl <- makeCluster(3)

  clusterEvalQ(cl, {library(moveLLR); library(ascr)})
  clusterExport(cl, varlist = c("disperse.start", "bincapt.session", "posts.full", "bearing.session", "animal.id.session", "survey.length"))

  maskSpecs.Results <- clusterApplyLB(cl, maskSpecs.List, \(x) {
    survey.buffer <- create.mask(posts.full, buffer = x$buffer)
    mask.dist.buffer <- lapply(1:length(survey.buffer), \(x) {ascr::distances(survey.buffer[[x]], posts.full[[x]])})
    expected.bearing.buffer <- lapply(1:length(survey.buffer), \(x) {t(ascr::bearings(posts.full[[x]], survey.buffer[[x]]))})
    pixel.size.buffer <- sapply(survey.buffer, attr, "area")

    buffer.model <- try(optim(disperse.start, fit_ascrDisperseSession, bin_capt = bincapt.session, mask = survey.buffer,
                              traps = posts.full, mask_dists = mask.dist.buffer, observed_bearings =  bearing.session,
                              expected_bearings = expected.bearing.buffer, pixel_area = pixel.size.buffer,
                              ID = animal.id.session, survey_length = rep(survey.length, length(bincapt.session)), threshold = 1e-4, method = "Nelder-Mead"))

    return (buffer.model)
  })

  stopCluster(cl)

  estimated.betas.disperse <- lapply(maskSpecs.Results, \(x) {x$par}) |>
    do.call(what = "rbind")

  colnames(estimated.betas.disperse) <- c("ell_D", "ell_lambda0", "ell_sigma", "logit_pb", "ell_sigmam", "ell_kappa")

  saveRDS(estimated.betas.disperse, "code/gibbons/03A-buffer-explore.rds")
} else {
  estimated.betas.disperse <- readRDS("code/gibbons/03A-buffer-explore.rds")
}

# -------------------------------------------------------------------------

maskSpecs <- cbind(maskSpecs, estimated.betas.disperse)
maskSpecs$D <- exp(maskSpecs$ell_D)
maskSpecs$lambda0 <- exp(maskSpecs$ell_lambda0)
maskSpecs$sigma <- exp(maskSpecs$ell_sigma)
maskSpecs$pb <- plogis(maskSpecs$logit_pb)
maskSpecs$sigmam <- exp(maskSpecs$ell_sigmam)
maskSpecs$kappa <- exp(maskSpecs$ell_kappa)

# -------------------------------------------------------------------------

layout(matrix(1:6, 2, 3))
plot(D ~ buffer, data = maskSpecs, type = "l", ylim = range(maskSpecs$D))
# lines(D ~ buffer, data = subset(maskSpecs, spacing == 200), type = "l", col = "red")
# lines(D ~ buffer, data = subset(maskSpecs, spacing == 250), type = "l", col = "blue")
# lines(D ~ buffer, data = subset(maskSpecs, spacing == 300), type = "l", col = "green")

plot(lambda0 ~ buffer, data = maskSpecs, type = "l", ylim = range(maskSpecs$lambda0))
# lines(lambda0 ~ buffer, data = subset(maskSpecs, spacing == 200), type = "l", col = "red")
# lines(lambda0 ~ buffer, data = subset(maskSpecs, spacing == 250), type = "l", col = "blue")
# lines(lambda0 ~ buffer, data = subset(maskSpecs, spacing == 300), type = "l", col = "green")

plot(sigma ~ buffer, data = maskSpecs, type = "l", ylim = range(maskSpecs$sigma))
# lines(sigma ~ buffer, data = subset(maskSpecs, spacing == 200), type = "l", col = "red")
# lines(sigma ~ buffer, data = subset(maskSpecs, spacing == 250), type = "l", col = "blue")
# lines(sigma ~ buffer, data = subset(maskSpecs, spacing == 300), type = "l", col = "green")

plot(pb ~ buffer, data = maskSpecs, type = "l", ylim = range(maskSpecs$pb))
# lines(pb ~ buffer, data = subset(maskSpecs, spacing == 200), type = "l", col = "red")
# lines(pb ~ buffer, data = subset(maskSpecs, spacing == 250), type = "l", col = "blue")
# lines(pb ~ buffer, data = subset(maskSpecs, spacing == 300), type = "l", col = "green")

plot(sigmam ~ buffer, data = maskSpecs, type = "l", ylim = range(maskSpecs$sigmam))
# lines(sigmam ~ buffer, data = subset(maskSpecs, spacing == 200), type = "l", col = "red")
# lines(sigmam ~ buffer, data = subset(maskSpecs, spacing == 250), type = "l", col = "blue")
# lines(sigmam ~ buffer, data = subset(maskSpecs, spacing == 300), type = "l", col = "green")

plot(kappa ~ buffer, data = maskSpecs, type = "l", ylim = range(maskSpecs$kappa))
# lines(kappa ~ buffer, data = subset(maskSpecs, spacing == 200), type = "l", col = "red")
# lines(kappa ~ buffer, data = subset(maskSpecs, spacing == 250), type = "l", col = "blue")
# lines(kappa ~ buffer, data = subset(maskSpecs, spacing == 300), type = "l", col = "green")

