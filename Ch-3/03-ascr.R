library(ascr)
source("./code/sonobuoys/00-aux-funcs.R")
load("./data/sonobuoys.RData")

surveyLength.session <- seaTrial.df %>%
  select(TrackedGroup, posixDate) %>%
  mutate(TrackedGroup = LETTERS[as.numeric(gsub("S([0-9]+)", "\\1", TrackedGroup))]) %>%
  filter(complete.cases(.)) %>%
  mutate(dateTime = lubridate::parse_date_time(posixDate, "dmYHM"),
         minute = lubridate::minute(dateTime)) %>%
  group_by(TrackedGroup) %>%
  summarise(duration = max(minute) - min(minute)) %>%
  mutate(duration.new = ifelse(duration == 0, 1, duration))

survey.length = sum(surveyLength.session$duration.new)

# Start Loc ---------------------------------------------------------------
start.hn <- fit.ascr(capt = collapsed.capt, traps = buoy.start.xy,
                     mask = mask.start, detfn = "hn", survey.length = survey.length)
start.hr <- fit.ascr(capt = collapsed.capt, traps = buoy.start.xy,
                     mask = mask.start, detfn = "hr", survey.length = survey.length)
start.th <- fit.ascr(capt = collapsed.capt, traps = buoy.start.xy,
                     mask = mask.start, detfn = "th", survey.length = survey.length)
start.hhn <- fit.ascr(capt = collapsed.capt, traps = buoy.start.xy,
                     mask = mask.start, detfn = "hhn", survey.length = survey.length)

# Start Abs Loc -----------------------------------------------------------
# start.abs.hn <- fit.ascr(capt = collapsed.capt, traps = buoy.abs.start.xy,
#                      mask = mask.abs.start, detfn = "hn")
# start.abs.hr <- fit.ascr(capt = collapsed.capt, traps = buoy.abs.start.xy,
#                      mask = mask.abs.start, detfn = "hr")
# start.abs.th <- fit.ascr(capt = collapsed.capt, traps = buoy.abs.start.xy,
#                      mask = mask.abs.start, detfn = "th")
# start.abs.hhn <- fit.ascr(capt = collapsed.capt, traps = buoy.abs.start.xy,
#                       mask = mask.abs.start, detfn = "hhn")

# Last Known Loc ----------------------------------------------------------
# last.hn <- fit.ascr(capt = collapsed.capt, traps = buoy.last.xy,
#                      mask = mask.last, detfn = "hn")
# last.hr <- fit.ascr(capt = collapsed.capt, traps = buoy.last.xy,
#                      mask = mask.last, detfn = "hr")
# last.th <- fit.ascr(capt = collapsed.capt, traps = buoy.last.xy,
#                      mask = mask.last, detfn = "th")
# last.hhn <- fit.ascr(capt = collapsed.capt, traps = buoy.last.xy,
#                       mask = mask.last, detfn = "hhn")

# Drift w/ Batches --------------------------------------------------------
survey.length.batch = surveyLength.session$duration.new

driftB.hn <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session,
                      mask = mask.obj.session, detfn = "hn", survey.length = survey.length.batch)
driftB.hr <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session,
                      mask = mask.obj.session, detfn = "hr", survey.length = survey.length.batch)
driftB.th <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session,
                      mask = mask.obj.session, detfn = "th", survey.length = survey.length.batch)
driftB.hhn <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session,
                       mask = mask.obj.session, detfn = "hhn", survey.length = survey.length.batch)

## hhn failed the converge ~ refitting with larger buffer using the same pixel area
mask.obj.session.mod <- seaTrial.df[, c("buoy.rel.x", "buoy.rel.y")] |>
  as.matrix() |>
  create.mask(buffer = 6000, spacing = sqrt(attr(mask.obj, "area") * 10000)) |>
  list() |>
  rep(length(capt.hist.session))
driftB.hhn <- fit.ascr(capt = capt.session, traps = buoy.ave.locs.session,
                       mask = mask.obj.session.mod, detfn = "hhn", survey.length = survey.length.batch)

# Drift w/out Batches -----------------------------------------------------
driftNB.hn <- fit.ascr(capt = capt.unbatched, traps = buoy.unbatched.xy,
                       mask = mask.unbatched, detfn = "hn")
driftNB.hr <- fit.ascr(capt = capt.unbatched, traps = buoy.unbatched.xy,
                       mask = mask.unbatched, detfn = "hr")
driftNB.th <- fit.ascr(capt = capt.unbatched, traps = buoy.unbatched.xy,
                       mask = mask.unbatched, detfn = "th")
driftNB.hhn <- fit.ascr(capt = capt.unbatched, traps = buoy.unbatched.xy,
                        mask = mask.unbatched, detfn = "hhn")

# Compiling for 07-inference.r --------------------------------------------

pattern.vec <- c("start\\.(hn|hr|th|hhn)", "start\\.abs\\.(hn|hr|th|hhn)", "driftB\\.(hn|hr|th|hhn)", "driftNB\\.(hn|hr|th|hhn)")
name.vec <- c("start", "start.abs", "drift.batch", "drift.unbatched")
names(name.vec) <- pattern.vec

ascr.models <- list()

for (whichPattern in pattern.vec) {
  ascr.models[[whichPattern]] <- lapply(ls(pattern = whichPattern), modelResPrinter, w = name.vec[whichPattern]) |>
    do.call(what = "rbind")
}

rbind(ascr.models[[1]], ascr.models[[2]], ascr.models[[3]]) |>
  saveRDS("./code/sonobuoys/ascr-models.RDS")
