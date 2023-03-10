source("./code/miss-v2/01-aux-funcs.R")

load(file = "code/miss-v2/lightfooti.RData")
lightfooti.admb = lightfooti
ftraps <- lightfooti.admb$traps

trap.dist = distances(ftraps, ftraps) |> mean()
trap.dist

survey.designs <- list(
  traps.line(2, trap.dist),
  traps.line(3, trap.dist),
  traps.line(4, trap.dist),
  traps.equilTriangle(trap.dist),
  traps.grid(4, trap.dist),
  traps.grid(9, trap.dist),
  traps.grid(25, trap.dist)
)

saveRDS(survey.designs, file = "./code/miss-v2/surveyDesigns.RDS")

simulation.reps <- 1000
set.seed(5353983)
cluster.seeds <- sample(1:1e7, simulation.reps, replace = FALSE)
write.csv(data.frame(seed = cluster.seeds, run.ids = rep(1:200, each = 5)), file = "./code/miss-v2/seeds.csv", row.names = FALSE)

# -------------------------------------------------------------------------

for (i in 1:length(survey.designs)) {
  it.filename <- paste0("./code/miss-v2/fascrNew_", i, ".csv")

  write.table(data.frame(
    theta = numeric(0), theta.name = character(0), survey = character(0),
    detfn = character(0), seed = numeric(0), aux = logical(0)),
    file = it.filename, row.names = FALSE, sep = ",")
}
