sigmaMove <- sqrt((30 * 10000) / pi) / 2.4477
nOccasion <- 3

full.df <- simul_data(
     density = 0.5, survey.region = aBounds, listening.post = listeningPosts,
     detfn.theta = c(lambda0 = 5, sigma = 750), kappa = 10,
     movement.theta = c(sigma_11 = sigmaMove, sigma_12 = 0, sigma_21 = 0, sigma_22 = sigmaMove)^2,
     call.Process = "Binomial", call.Process.Theta = c(nOccasion, 2/3)
 )


# -------------------------------------------------------------------------

mask <- create.mask(listeningPosts, buffer = 4500)

plot(ac.y ~ ac.x, data = full.df, xlim = c(-10000, 10000), col = "grey50")
points(mask, col = "grey90")

points(listeningPosts, pch  = 16 , col = "red")
points(y ~ x, data = full.df %>% mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>% filter(detected > 0), pch = 16, col = "grey50")
points(ac.y ~ ac.x, data = full.df %>% mutate(detected = rowSums(.[grep("detect\\.", names(.))])) %>% filter(detected > 0), pch = 16)

# -------------------------------------------------------------------------


mvnorm_nll <- function(theta, y, ac) {
  sum(-log(mvtnorm::dmvnorm(x = y, mean = ac, sigma = diag(exp(theta), 2)^2)))
}

sigma.hat <- unique(full.df$groupID) |> length() |> numeric()

for (id in unique(full.df$groupID)) {
  group.id <- subset(full.df, groupID == id)

  sigma.est <- optim(log(150), mvnorm_nll, y = as.matrix(group.id[, c("x", "y")]), ac = group.id[, c("ac.x", "ac.y")] |> unique() |> as.numeric(),
                     method = "L-BFGS-B")

  sigma.est$par |> exp() -> sigma.hat[id]
}

terri.hat <- (2.4477 * sigma.hat)^2 * pi / 10000

hist(terri.hat)

summary(terri.hat)

