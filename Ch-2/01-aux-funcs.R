sim.capt.wrap <- function(whichScenario, traps) {
  mask <- switch(whichScenario,
                     oascr = read.traps(data = as.data.frame(traps)) |> secr::make.mask(buffer = 200) |> as.matrix(),
                     fascr = ascr::create.mask(traps, buffer = 40),
                     gascr = ascr::create.mask(traps, buffer = 6000))

  list(
    ascr::sim.capt(traps = traps, mask = mask, infotypes = sim.theta[[whichScenario]]$infotypes,
           pars = sim.theta[[whichScenario]]$theta, detfn = sim.theta[[whichScenario]]$detfn,
           ss.opts = sim.theta[[whichScenario]]$ss.opts),
    mask
  )
}

# -------------------------------------------------------------------------

traps.line <- function(n.traps, spacing, horizontal = TRUE) {
  if (horizontal) {
    return (cbind(x = seq(0, spacing * (n.traps - 1), length.out = n.traps), y = rep(0, n.traps)))
  }

  return (cbind(x = rep(0, n.traps), y = seq(0, spacing * (n.traps - 1), length.out = n.traps)))
}

traps.grid <- function(n.traps, spacing) {
  if (sqrt(n.traps) != round(sqrt(n.traps))) {
    stop("Not a valid grid")
  }

  n.nodes <- sqrt(n.traps)

  expand.grid(x = seq(0, spacing * (n.nodes - 1), length.out = n.nodes),
              y = seq(0, spacing * (n.nodes - 1), length.out = n.nodes)) |>
    as.matrix()
}

traps.equilTriangle <- function(spacing) {
  cbind(x = seq(0, spacing, length.out = 3), y = c(0, sqrt(3) / 2 * spacing, 0))
}
