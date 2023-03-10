library(readxl); library(ascr); library(secr); library(Rcpp); library(RcppArmadillo)

load("code/gibbons/01-data-prep.RData")
sourceCpp("code/gibbons/cpp-Files/ascrStack.cpp"); sourceCpp("code/gibbons/cpp-Files/ascrDisperse.cpp")

# -------------------------------------------------------------------------

benchmark.fit <- fit.ascr(capt.bearings, traps = posts.full, detfn = "hhn", mask = survey.mask)

# -------------------------------------------------------------------------

stack.start <- c(log(0.005), log(2), log(650), qlogis(0.5), log(8.28))
stack.model <- optim(stack.start, fit_ascrStackSession, bin_capt = bincapt.session, mask = survey.session,
                     traps = posts.full, mask_dists = mask.dist.session, observed_bearings =  bearing.session,
                     expected_bearings = expected.bearing.session, pixel_area = pixel.size.session,
                     ID = animal.id.session, survey_length = rep(survey.length, length(bincapt.session)), bernoulli = TRUE,
                     hessian = TRUE, method = "L-BFGS-B")

# -------------------------------------------------------------------------

disperse.start <- c(log(0.005), log(2), log(650), qlogis(0.5), log(550), log(8.28))
disperse.model <- optim(disperse.start, fit_ascrDisperseSession, bin_capt = bincapt.session, mask = survey.session,
                        traps = posts.full, mask_dists = mask.dist.session, observed_bearings =  bearing.session,
                        expected_bearings = expected.bearing.session, pixel_area = pixel.size.session,
                        ID = animal.id.session, survey_length = rep(survey.length, length(bincapt.session)), threshold = 1e-4,
                        trace = TRUE, hessian = TRUE, method = "L-BFGS-B")

# -------------------------------------------------------------------------

save(benchmark.fit, stack.model, disperse.model, file = "code/gibbons/03-fit-models.RData")
