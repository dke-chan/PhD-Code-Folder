full.set <- expand.grid(detfn = c("hn", "hn.f", "hr", "th", "hhn"),
                        covariate.1 = c("coast", NA),
                        covariate.2 = c("depth", NA),
                        k.1 = 1:4,
                        k.2 = 1:4,
                        stringsAsFactors = FALSE) |>
  subset(!is.na(covariate.1) | !is.na(covariate.2)) |>
  {\(x){
    y <- x[rep(1:nrow(x), 2), ]
    y$basis <- rep(c("tp", "cr"), each = nrow(x))
    y$covariate.input <- ifelse(
      !is.na(y$covariate.1) & !is.na(y$covariate.2),
      "both",
      ifelse(!is.na(y$covariate.1), y$covariate.1, y$covariate.2)
    )
    z <- y[-which(y$covariate.input != "both" & y$k.1 != y$k.2), ]

    return (z)
  }}() |>
  rbind(data.frame(
    detfn = c("hn", "hn.f", "hr", "th", "hhn"),
    covariate.1 = NA,
    covariate.2 = NA,
    k.1 = 0,
    k.2 = 0,
    basis = NA,
    covariate.input = "none"
  ))

write.csv(full.set, "./code/whales/scr-to-fit.csv", row.names = FALSE)
