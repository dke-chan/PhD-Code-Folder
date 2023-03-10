# library(spatstat.geom, include.only = c("pairdist", "crossdist"))
# library(spatstat.utils, include.only = "matrowany")
# library(CircStats, include.only = "rvm")
# library(mvtnorm, include.only = "rmvnorm")
# library(dplyr, include.only = c("filter", "between"))
# library(ascr, include.only = "bearings")

simul_data <- function (
  density = 0.75,
  survey.region = expand.grid(x = seq(-5000, 5000, by = 250), y = seq(-5000, 5000, by = 250)),
  listening.post = data.frame(x = seq(-1000, 1000, by = 1000), y = 0),
  occasions = 1,
  detfn = "hhn", detfn.theta = c(lambda0 = 5, sigma = 1250),
  kappa = 36.22,
  movement.dist = "mvn", movement.theta = c(sigma_11 = 1000, sigma_12 = 0, sigma_21 = 0, sigma_22 = 1000)^2,
  call.Process = "Poisson", call.Process.Theta = 1.25,
  point.Process = "Poisson", maternII.bound = 500) {

  ## Note that everything is either in terms of metres or converted to be in terms of metres
  ##
  ## density: Expected number of animals/groups per km^2
  ## occasions: The number of sampling occasions
  ## survey.region: The survey region
  ## listening.post: The coordinates of the detectors in the survey region
  ## kappa: The von Mises' scale parameter for modelling the bearing error
  ## detfn: Detection function ~ See ?secr::detectfn for more information
  ## detfn.theta: The parameters of the detection function
  ## movement.dist: The distribution for modelling the movement about territory centres
  ## movement.theta: The parameters of the movement distribution
  ## call.Process: The distribution which models the number of calls made per survey
  ## call.Process.Theta: The parameters of the number of calls distribution
  ## point.Process: The underlying point process for the true distribution of activity centres
  ## maternII.bound: The boundary for determining whether a territory centre should be removed

  # occasions not actually functioning as intended btw for Binomial / Bernoulli scenarios

  density.m2 <- density / 1e+6 ## Convert to animals/groups per m^2

  ## Generate territory centres
  if (point.Process == "Poisson") {
    area <- diff(range(survey.region$x)) * diff(range(survey.region$y)) ## area of a square

    true.N <- rpois(1, density.m2 * area) ## the true number of groups

    true.ActivityCentre <- data.frame(x = runif(true.N, min(survey.region$x), max(survey.region$x)), ## the true centroid for each group
                                      y = runif(true.N, min(survey.region$y), max(survey.region$y)))
  } else if (point.Process == "MaternII") {
    buffer <- (diff(range(survey.region$x)) + diff(range(survey.region$y))) / 2 ## buffer range in m
    buff.survey.region <- expand.grid(x = seq(min(survey.region$x) - buffer, max(survey.region$x) + buffer, by = 250), y = seq(min(survey.region$y) - buffer, max(survey.region$y) + buffer, by = 250))

    buff.area <- diff(range(buff.survey.region$x)) * diff(range(buff.survey.region$y)) ## area of a square

    true.N.Star <- rpois(1, density.m2 * buff.area) ## the true number of groups

    ## Adapated from spatstat.core::rMaternInhibition
    true.ActivityCentre.Star <- data.frame(x = runif(true.N.Star, min(buff.survey.region$x), max(buff.survey.region$x)), ## the true centroid for each group
                                           y = runif(true.N.Star, min(buff.survey.region$y), max(buff.survey.region$y)),
                                           age = sample(seq_len(true.N.Star), true.N.Star, replace = FALSE))

    bound.Met <- (spatstat.geom::pairdist(true.ActivityCentre.Star[, c("x", "y")])^2) <= (maternII.bound^2)
    earlier.Age <- outer(true.ActivityCentre.Star$age, true.ActivityCentre.Star$age, ">")
    conflict.list <- bound.Met & earlier.Age
    toRemove <- spatstat.utils::matrowany(conflict.list)

    true.ActivityCentre <- true.ActivityCentre.Star[!toRemove, c("x", "y")]
    true.ActivityCentre <- dplyr::filter(true.ActivityCentre,
                                         dplyr::between(x, min(survey.region$x), max(survey.region$x)),
                                         dplyr::between(y, min(survey.region$y), max(survey.region$y)))

    true.N <- nrow(true.ActivityCentre) ## the true number of groups
  }

  ## Set-up the data.frame to be returned to the end user
  no.of.detectors <- nrow(listening.post)
  simul.frame.colnames <- c("x", "y", "ac.x", "ac.y", "groupID", "callNum", "occasionID")

  for (variable in c("dist", "bearing", "detect", "ob.bearing")) {
    simul.frame.colnames <- c(simul.frame.colnames, paste(variable, 1:no.of.detectors, sep = "."))
  }

  simul.frame.empty <- vector("list", length(simul.frame.colnames))
  simul.frame.empty <- lapply(1:length(simul.frame.colnames), function (x) 0)

  simul.frame <- do.call(data.frame, simul.frame.empty)
  colnames(simul.frame) <- simul.frame.colnames
  simul.frame <- simul.frame[-1, ]

  ## Detection function definitions...
  hhn.fn <- function(d, theta) {
    lambda0 <- theta[1]
    sigma.2 <- theta[2]^2

    lambda.d <- lambda0 * exp(-d^2 / (2 * sigma.2))
    1 - exp(-lambda.d)
  }

  ## Auxiliary calls...
  movement.theta.mat <- matrix(movement.theta, ncol = 2) ## "Typecasting"

  ## For each occasion...
  for (occasion.id in 1:occasions) {
    ## The number of calls process
    call.Process.Vector <- numeric(true.N)
    if (call.Process == "Poisson") {
      call.Process.Vector <- rpois(true.N, call.Process.Theta)
    } else if (call.Process == "Binomial") {
      call.Process.Vector <- rbinom(true.N, call.Process.Theta[1], call.Process.Theta[2])
    } else if (call.Process == "Bernoulli") {
      call.Process.Vector <- sapply(1:true.N, function (x) {
        sample(0:1, size = call.Process.Theta[1], replace = TRUE, prob = c(1 - call.Process.Theta[2], call.Process.Theta[2])) |> sum()
      })
    }

    ## Populate the data.frame!
    for (id in 1:true.N) {
      append.Frame <- true.ActivityCentre[rep(id, max(call.Process.Vector[id], 1)), ] + mvtnorm::rmvnorm(max(call.Process.Vector[id], 1), sigma = movement.theta.mat)

      if ((call.Process == "Binomial" | call.Process == "Bernoulli") & call.Process.Vector[id] < call.Process.Theta[1]) {
        append.Frame <- rbind(
          append.Frame,
          true.ActivityCentre[rep(id, call.Process.Theta[1] - nrow(append.Frame)), ] + mvtnorm::rmvnorm(call.Process.Theta[1] - nrow(append.Frame), sigma = movement.theta.mat)
        )
      }

      append.Frame[, paste0("ac.", c("x", "y"))] <- true.ActivityCentre[id, ]
      append.Frame[, "groupID"] <- id
      append.Frame[, "callNum"] <- call.Process.Vector[id]
      append.Frame[, "occasionID"] <- occasion.id

      append.Frame[, paste0("dist.", 1:no.of.detectors)] <- spatstat.geom::crossdist(append.Frame$x, append.Frame$y, listening.post$x, listening.post$y)
      append.Frame[, paste0("bearing.", 1:no.of.detectors)] <- t(ascr::bearings(as.matrix(listening.post), as.matrix(append.Frame[, c('x', 'y')])))

      if (call.Process.Vector[id] > 0) {
        append.Frame[1:call.Process.Vector[id], paste0("detect.", 1:no.of.detectors)] <- rbinom(call.Process.Vector[id] * no.of.detectors, 1, hhn.fn(unlist(append.Frame[1:call.Process.Vector[id], paste0("dist.", 1:no.of.detectors)], use.names = FALSE), detfn.theta))
        append.Frame[1:call.Process.Vector[id], paste0("ob.bearing.", 1:no.of.detectors)] <- sapply(unlist(append.Frame[1:call.Process.Vector[id], paste0("bearing.", 1:no.of.detectors)], use.names = FALSE), CircStats::rvm, n = 1, k = kappa)
      } else {
        append.Frame[, paste0("detect.", 1:no.of.detectors)] <- 0
        append.Frame[, paste0("ob.bearing.", 1:no.of.detectors)] <- 0
      }

      simul.frame <- rbind(simul.frame, append.Frame)
    }
  }


  ## For each occasion...
  if ((call.Process == "Binomial" | call.Process == "Bernoulli") & call.Process.Theta[1] > 1) {
    for (id in unique(simul.frame$groupID)) {
      simul.frame$occasionID[which(simul.frame$groupID == id)] <- sample(1:call.Process.Theta[1])
    }
  }

  return (simul.frame)
}

simul_session <- function (
  n.session = 2, listening.post.session = rep(list(data.frame(x = seq(-1000, 1000, by = 1000), y = 0)), n.session),
  ...
) {
  ## Use a dummy post for simul_data
  arg.list <- list(...)
  arg.list$listening.post <- data.frame(x = 0, y = 0)

  ## Load in the default values from simul_data if none are provided
  simul_data.defaults <- formals(simul_data)
  simul_data.arguments <- names(simul_data.defaults)

  for (label in simul_data.arguments[which(simul_data.arguments != "listening.post")]) {
    if (is.null(arg.list[[label]])) {
      if (label == "detfn.theta") {
        arg.list[[label]] <- c(lambda0 = 5, sigma = 1250)
      } else if (label == "movement.theta") {
        arg.list[[label]] <- c(sigma_11 = 1000, sigma_12 = 0, sigma_21 = 0, sigma_22 = 1000)^2
      } else {
        arg.list[[label]] <- simul_data.defaults[[label]]
      }
    }
  }

  ## Simulate the movement of animals with simul_data()
  all.animals <- do.call(simul_data, arg.list)
  all.animals.callIndex <- !is.na(all.animals$detect.1)

  ## Detection function definitions...
  hhn.fn <- function(d, theta) {
    lambda0 <- theta[1]
    sigma.2 <- theta[2]^2

    lambda.d <- lambda0 * exp(-d^2 / (2 * sigma.2))
    1 - exp(-lambda.d)
  }

  ## Simulate the capture of each animal wrt to each session!
  toReturn <- rep(list(subset(all.animals, select = c("x", "y", "ac.x", "ac.y", "groupID", "callNum", "occasionID"))), n.session)

  for (sessionIndex in 1:n.session) {
    no.of.detectors <- nrow(listening.post.session[[sessionIndex]])
    toReturn[[sessionIndex]][, paste0("dist.", 1:no.of.detectors)] <- spatstat.geom::crossdist(toReturn[[sessionIndex]]$x, toReturn[[sessionIndex]]$y, listening.post.session[[sessionIndex]]$x, listening.post.session[[sessionIndex]]$y)
    toReturn[[sessionIndex]][, paste0("bearing.", 1:no.of.detectors)] <- t(ascr::bearings(as.matrix(listening.post.session[[sessionIndex]]), as.matrix(toReturn[[sessionIndex]][, c('x', 'y')])))

    toReturn[[sessionIndex]][, paste0("detect.", 1:no.of.detectors)] <- rbinom(length(all.animals.callIndex) * no.of.detectors, 1, hhn.fn(unlist(toReturn[[sessionIndex]][, paste0("dist.", 1:no.of.detectors)], use.names = FALSE), arg.list[["detfn.theta"]]))
    toReturn[[sessionIndex]][!all.animals.callIndex, paste0("detect.", 1:no.of.detectors)] <- NA
    toReturn[[sessionIndex]][, paste0("ob.bearing.", 1:no.of.detectors)] <- sapply(unlist(toReturn[[sessionIndex]][, paste0("bearing.", 1:no.of.detectors)], use.names = FALSE), CircStats::rvm, n = 1, k = arg.list[["kappa"]])
  }

  return (toReturn)
}

simul_pixel <- function (
  iter = 100, x = 0, y = 0,
  listening.post = data.frame(x = seq(-1000, 1000, by = 1000), y = 0),
  detfn = "hhn", detfn.theta = c(lambda0 = 5, sigma = 1250),
  kappa = 36.22,
  occasions = 1,
  movement.dist = "mvn", movement.theta = c(sigma_11 = 1000, sigma_12 = 0, sigma_21 = 0, sigma_22 = 1000)^2,
  call.Process = "Poisson", call.Process.Theta = 1.25,
  point.Process = "Poisson", maternII.bound = 500) {

  ## Set-up the data.frame to be returned to the end user
  no.of.detectors <- nrow(listening.post)
  simul.frame.colnames <- c("x", "y", "ac.x", "ac.y", "groupID", "callNum", "occasionID")

  for (variable in c("dist", "bearing", "detect", "ob.bearing")) {
    simul.frame.colnames <- c(simul.frame.colnames, paste(variable, 1:no.of.detectors, sep = "."))
  }

  simul.frame.empty <- vector("list", length(simul.frame.colnames))
  simul.frame.empty <- lapply(1:length(simul.frame.colnames), function (x) 0)

  simul.frame <- do.call(data.frame, simul.frame.empty)
  colnames(simul.frame) <- simul.frame.colnames
  simul.frame <- simul.frame[-1, ]

  ## Detection function definitions...
  hhn.fn <- function(d, theta) {
    lambda0 <- theta[1]
    sigma.2 <- theta[2]^2

    lambda.d <- lambda0 * exp(-d^2 / (2 * sigma.2))
    1 - exp(-lambda.d)
  }

  ##
  true.ActivityCentre <- data.frame(x = rep(x, iter), ## the true centroid for each group
                                    y = rep(y, iter))

  ## Auxiliary calls...
  movement.theta.mat <- matrix(movement.theta, ncol = 2) ## "Typecasting"

  ## For each occasion...
  for (occasion.id in 1:occasions) {
    ## The number of calls process
    call.Process.Vector <- numeric(iter)
    if (call.Process == "Poisson") {
      call.Process.Vector <- rpois(iter, call.Process.Theta)
    } else if (call.Process == "Binomial") {
      call.Process.Vector <- rbinom(iter, call.Process.Theta[1], call.Process.Theta[2])
    } else if (call.Process == "Bernoulli") {
      call.Process.Vector <- sapply(1:iter, function (x) {
        sample(0:1, size = call.Process.Theta[1], replace = TRUE, prob = c(1 - call.Process.Theta[2], call.Process.Theta[2])) |> sum()
      })
    }

    ## Populate the data.frame!
    for (id in 1:iter) {
      append.Frame <- true.ActivityCentre[rep(id, max(call.Process.Vector[id], 1)), ] + mvtnorm::rmvnorm(max(call.Process.Vector[id], 1), sigma = movement.theta.mat)

      if ((call.Process == "Binomial" | call.Process == "Bernoulli") & call.Process.Vector[id] < call.Process.Theta[1]) {
        append.Frame <- rbind(
          append.Frame,
          true.ActivityCentre[rep(id, call.Process.Theta[1] - nrow(append.Frame)), ] + mvtnorm::rmvnorm(call.Process.Theta[1] - nrow(append.Frame), sigma = movement.theta.mat)
        )
      }

      append.Frame[, paste0("ac.", c("x", "y"))] <- true.ActivityCentre[id, ]
      append.Frame[, "groupID"] <- id
      append.Frame[, "callNum"] <- call.Process.Vector[id]
      append.Frame[, "occasionID"] <- occasion.id

      append.Frame[, paste0("dist.", 1:no.of.detectors)] <- spatstat.geom::crossdist(append.Frame$x, append.Frame$y, listening.post$x, listening.post$y)
      append.Frame[, paste0("bearing.", 1:no.of.detectors)] <- t(ascr::bearings(as.matrix(listening.post), as.matrix(append.Frame[, c('x', 'y')])))

      if (call.Process.Vector[id] > 0) {
        append.Frame[1:call.Process.Vector[id], paste0("detect.", 1:no.of.detectors)] <- rbinom(call.Process.Vector[id] * no.of.detectors, 1, hhn.fn(unlist(append.Frame[1:call.Process.Vector[id], paste0("dist.", 1:no.of.detectors)], use.names = FALSE), detfn.theta))
        append.Frame[1:call.Process.Vector[id], paste0("ob.bearing.", 1:no.of.detectors)] <- sapply(unlist(append.Frame[1:call.Process.Vector[id], paste0("bearing.", 1:no.of.detectors)], use.names = FALSE), CircStats::rvm, n = 1, k = kappa)
      } else {
        append.Frame[, paste0("detect.", 1:no.of.detectors)] <- 0
        append.Frame[, paste0("ob.bearing.", 1:no.of.detectors)] <- 0
      }

      simul.frame <- rbind(simul.frame, append.Frame)
    }
  }

  ## For each occasion...
  if ((call.Process == "Binomial" | call.Process == "Bernoulli") & call.Process.Theta[1] > 1) {
    for (id in unique(simul.frame$groupID)) {
      simul.frame$occasionID[which(simul.frame$groupID == id)] <- sample(1:call.Process.Theta[1])
    }
  }

  return (simul.frame)
}
