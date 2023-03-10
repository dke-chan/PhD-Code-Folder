library(ascr)
source("./code/sonobuoys/00-aux-funcs.R")
load("./data/sonobuoys.RData")

# Start Loc ---------------------------------------------------------------
# start.hn <- fit.ascr(capt = collapsed.mrds, traps = buoy.start.xy,
#                      mask = mask.start, detfn = "hn")
# start.hr <- fit.ascr(capt = collapsed.mrds, traps = buoy.start.xy,
#                      mask = mask.start, detfn = "hr")
# start.th <- fit.ascr(capt = collapsed.mrds, traps = buoy.start.xy,
#                      mask = mask.start, detfn = "th")
# start.hhn <- fit.ascr(capt = collapsed.mrds, traps = buoy.start.xy,
#                      mask = mask.start, detfn = "hhn")

# Last Known Loc ----------------------------------------------------------
# last.hn <- fit.ascr(capt = collapsed.mrds, traps = buoy.last.xy,
#                      mask = mask.last, detfn = "hn")
# last.hr <- fit.ascr(capt = collapsed.mrds, traps = buoy.last.xy,
#                      mask = mask.last, detfn = "hr")
# last.th <- fit.ascr(capt = collapsed.mrds, traps = buoy.last.xy,
#                      mask = mask.last, detfn = "th")
# last.hhn <- fit.ascr(capt = collapsed.mrds, traps = buoy.last.xy,
#                       mask = mask.last, detfn = "hhn")

# Drift w/ Batches --------------------------------------------------------
driftB.hn <- fit.ascr(capt = capt.session.mrds, traps = buoy.ave.locs.session,
                      mask = mask.obj.session, detfn = "hn")
driftB.hr <- fit.ascr(capt = capt.session.mrds, traps = buoy.ave.locs.session,
                      mask = mask.obj.session, detfn = "hr")
driftB.th <- fit.ascr(capt = capt.session.mrds, traps = buoy.ave.locs.session,
                      mask = mask.obj.session, detfn = "th")
driftB.hhn <- fit.ascr(capt = capt.session.mrds, traps = buoy.ave.locs.session,
                       mask = mask.obj.session, detfn = "hhn")

# hr failed the converge ~ refitting with larger buffer using the same pixel area
mask.obj.session.mod <- seaTrial.df[, c("buoy.rel.x", "buoy.rel.y")] |>
  as.matrix() |>
  create.mask(buffer = 9000, spacing = sqrt(attr(mask.obj, "area") * 10000)) |>
  list() |>
  rep(length(capt.hist.session))
driftB.hr <- fit.ascr(capt = capt.session.mrds, traps = buoy.ave.locs.session,
                      mask = mask.obj.session.mod, detfn = "hr")

# Drift w/out Batches -----------------------------------------------------
driftNB.hn <- fit.ascr(capt = capt.unbatched.mrds, traps = buoy.unbatched.xy,
                       mask = mask.unbatched, detfn = "hn")
driftNB.hr <- fit.ascr(capt = capt.unbatched.mrds, traps = buoy.unbatched.xy,
                       mask = mask.unbatched, detfn = "hr")
driftNB.th <- fit.ascr(capt = capt.unbatched.mrds, traps = buoy.unbatched.xy,
                       mask = mask.unbatched, detfn = "th")
driftNB.hhn <- fit.ascr(capt = capt.unbatched.mrds, traps = buoy.unbatched.xy,
                        mask = mask.unbatched, detfn = "hhn")

# Compiling for 07-inference.r --------------------------------------------

mrds.models <- lapply(ls(pattern = "driftNB\\.(hn|hr|th|hhn)"), modelResPrinter, w = "drift.unbatched") |>
  do.call(what = "rbind")

saveRDS(mrds.models, "./code/sonobuoys/mrds-models.RDS")
