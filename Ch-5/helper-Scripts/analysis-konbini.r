thetaExtract <- function (data, whichModel, whichDesign, truth, exclude.outliers = TRUE, include.derived = FALSE) {
  dataFil <- data[which(data$model == whichModel & data$design == whichDesign), ]

  if (unique(dataFil$callProcess) == "Binomial") {
    coefNames <- switch(whichModel,
                        ascr = c("D_c", "lambda0", "sigma", "kappa"),
                        ascrStack = c("ell_D", "ell_lambda0", "ell_sigma", "logit_pb", "ell_kappa"),
                        ascrDisperse = c("ell_D", "ell_lambda0", "ell_sigma", "logit_pb", "ell_sigmam", "ell_kappa")
    )

    # backTransform <- switch(whichModel, ascr = FALSE, ascrStack = TRUE, ascrDisperse = TRUE)
  }

  nameMap_tX <- \(x) {
    switch(x,
           D_c = "ell_D_c",
           lambda0 = "ell_lambda0",
           sigma = "ell_sigma",
           kappa = "ell_kappa",
           x)
  }

  nameMap_btX <- \(x) {
    switch(x,
           ell_D = "D",
           ell_lambda0 = "lambda0",
           ell_sigma = "sigma",
           logit_pb = "pb",
           ell_sigmam = "sigmam",
           ell_kappa = "kappa",
           x)
  }

  transform_X <- \(thetaName, x) {
    switch(thetaName,
           D_c = log(x),
           lambda0 = log(x),
           sigma = log(x),
           kappa = log(x),
           x)
  }

  backTransform_X <- \(thetaName, x) {
    switch(thetaName,
           ell_D = exp(x),
           ell_lambda0 = exp(x),
           ell_sigma = exp(x),
           logit_pb = plogis(x),
           ell_sigmam = exp(x),
           ell_kappa = exp(x),
           D_c = exp(x), ## ...
           lambda0 = exp(x),
           sigma = exp(x),
           kappa = exp(x),
           x)
  }

  backTransform_SE <- \(thetaName, x) {
    switch(thetaName,
           ell_D = exp(x),
           ell_lambda0 = exp(x),
           ell_sigma = exp(x),
           logit_pb = dlogis(x),
           ell_sigmam = exp(x),
           ell_kappa = exp(x),
           x)
  }

  toReturn <- lapply(coefNames, \(x) {
    y <- dataFil[which(dataFil$coefName == x), ]
    whichMappedName <- nameMap_btX(x)

    toChild <- list()

    toChild$estimates <- transform_X(x, y$theta)
    toChild$se <- y$se ## Returns linked par std.errs

    # if (backTransform) {
      # toChild$estimates <- backTransform_X(x, toChild$estimates)
      # toChild$se <- sqrt(backTransform_SE(x, y$theta)^2 * toChild$se^2)
    # }

    if (exclude.outliers) {
      if (whichModel != "ascr") {
        a <- backTransform_X(x, toChild$estimates)
      } else {
        a <- toChild$estimates
      }

      toChild.Outliers <- boxplot.stats(a)$out
      Outliers.indices <- a %in% toChild.Outliers
      toChild$estimates[Outliers.indices] <- NA
      toChild$se[Outliers.indices] <- NA
    }

    toChild$lwrWalds <- toChild$estimates - qnorm(0.975) * toChild$se
    toChild$uprWalds <- toChild$estimates + qnorm(0.975) * toChild$se

    toChild$estimates <- backTransform_X(x, toChild$estimates)
    toChild$lwrWalds <- backTransform_X(x, toChild$lwrWalds)
    toChild$uprWalds <- backTransform_X(x, toChild$uprWalds)

    toChild$mean <- mean(toChild$estimates, na.rm = TRUE) |> unname()
    toChild$median <- median(toChild$estimates, na.rm = TRUE) |> unname()

    toChild$capture <- 100 * mean(toChild$lwrWalds <= truth[whichMappedName] & truth[whichMappedName] <= toChild$uprWalds, na.rm = TRUE) |> unname()
    toChild$bias <- mean(toChild$estimates, na.rm = TRUE) - truth[whichMappedName] |> unname()
    toChild$percentBias <- 100 * toChild$bias / truth[whichMappedName] |> unname()
    toChild$cv <- sd(toChild$estimates, na.rm = TRUE) / mean(toChild$estimates, na.rm = TRUE) * 100 |> unname()

    return (toChild)
  })

  toReturn <- setNames(toReturn, sapply(coefNames, nameMap_btX))

  if (include.derived & (whichModel == "ascrStack" | whichModel == "ascrDisperse")) {
    toReturn[["D_c"]] <- list()

    toReturn[["D_c"]]$estimates <- toReturn[["D"]]$estimates * toReturn[["pb"]]$estimates * truth["nc"]

    D_c_se <- \(D.theta, pb.theta, D.se, pb.se, cov.D.pb, nc) {
      toReturn <- numeric(length(D.theta))
      # print(length(D.theta))
      for (i in 1:length(D.theta)) {
        vcov.hat <- matrix(c(D.se[i]^2, cov.D.pb[i], cov.D.pb[i], pb.se[i]^2), 2, 2)

        G.1 <- nc * c(D.theta[i] * pb.theta[i], D.theta[i] * exp(qlogis(pb.theta[i])) / (exp(qlogis(pb.theta[i])) + 1)^2) |> matrix(1, 2)

        toReturn[i] <- sqrt(G.1 %*% vcov.hat %*% t(G.1))
      }

      return (toReturn)
    }

    toReturn[["D_c"]]$se <- D_c_se(toReturn[["D"]]$estimates, toReturn[["pb"]]$estimates, toReturn[["D"]]$se, toReturn[["pb"]]$se,
                                   dataFil$se[which(dataFil$coefName == "cov_ell_D_logit_pb")],
                                   truth["nc"])

    toReturn[["D_c"]]$lwrWalds <- toReturn[["D_c"]]$estimates - qnorm(0.975) * toReturn[["D_c"]]$se
    toReturn[["D_c"]]$uprWalds <- toReturn[["D_c"]]$estimates + qnorm(0.975) * toReturn[["D_c"]]$se

    toReturn[["D_c"]]$mean <- mean(toReturn[["D_c"]]$estimates, na.rm = TRUE) |> unname()
    toReturn[["D_c"]]$median <- median(toReturn[["D_c"]]$estimates, na.rm = TRUE) |> unname()

    toReturn[["D_c"]]$capture <- 100 * mean(toReturn[["D_c"]]$lwrWalds <= truth["nc"] * truth["D"] * truth["pb"] & truth["nc"] * truth["D"] * truth["pb"] <= toReturn[["D_c"]]$uprWalds, na.rm = TRUE) |> unname()
    toReturn[["D_c"]]$bias <- mean(toReturn[["D_c"]]$estimates, na.rm = TRUE) - truth["D_c"] |> unname()
    toReturn[["D_c"]]$percentBias <- 100 * toReturn[["D_c"]]$bias / truth["D_c"] |> unname()
    toReturn[["D_c"]]$cv <- sd(toReturn[["D_c"]]$estimates, na.rm = TRUE) / mean(toReturn[["D_c"]]$estimates, na.rm = TRUE) * 100 |> unname()
  }

  if (include.derived & whichModel == "ascrDisperse") {

  }


  return (toReturn)
}


# -------------------------------------------------------------------------

thetaBoxplot <- function (data, whichParameter, whichModels = c("ascr", "ascrStack", "ascrDisperse"), whichDesigns, truths, exclude.outliers = TRUE,
                          modelCol = NULL, ...) {
  whichParameter <- whichParameter[1]

  refGrid <- expand.grid(model = whichModels, design = whichDesigns, stringsAsFactors = FALSE)
  if (is.null(modelCol)) {
    refGrid$Col <- setNames(RColorBrewer::brewer.pal(3, "Dark2"), c("ascr", "ascrStack", "ascrDisperse"))[refGrid$model]
  } else {

    refGrid$Col <- modelCol
  }

  processed <- split(refGrid, ~ model + design) |>
    lapply(\(x) {
      y <- thetaExtract(data, x$model, x$design, truths[[x$design]], exclude.outliers, include.derived = TRUE)

      return (y[[whichParameter]]$estimates)
    })

  maxLength <- sapply(processed, length) |>
    max()

  work.DF <- rep(NA, times = maxLength * length(processed)) |>
    matrix(nrow = maxLength, ncol = length(processed)) |>
    as.data.frame()

  for (i in 1:length(processed)) {
    work.DF[1:length(processed[[i]]), i] <- processed[[i]]
  }

  boxplot(work.DF, col = refGrid$Col, xlab = "", ylab = "", names = rep("", nrow(refGrid)), axes = FALSE,
          outpch = 21, outbg = refGrid$Col, ...)
}

# -------------------------------------------------------------------------

thetaResults <- function (data, whichModels = c("ascr", "ascrStack", "ascrDisperse"), whichDesigns, truths, exclude.outliers = TRUE,
                          clean = TRUE, derived = TRUE) {
  processed <- expand.grid(model = whichModels, design = whichDesigns) |>
    apply(1, \(x) {
      y <- thetaExtract(data, x["model"], x["design"], truths[[x["design"]]], exclude.outliers, include.derived = derived)

      z <- data.frame(Model = c(x["model"], rep(ifelse(clean, "", x["model"]), length(y) - 1)),
                      Design = c(x["design"], rep(ifelse(clean, "", x["design"]), length(y) - 1)),
                      Theta = names(y),
                      Truth = truths[[x["design"]]][names(y)],
                      Mean = sapply(y, \(a) a$mean),
                      Median = sapply(y, \(a) a$median),
                      Bias = sapply(y, \(a) a$bias),
                      pctBias = sapply(y, \(a) a$percentBias),
                      CV = sapply(y, \(a) a$cv),
                      waldCapture = sapply(y, \(a) a$capture)
      )

      rownames(z) <- NULL

      return (z)
    }) |> do.call(what = "rbind")

  return (processed)
}
