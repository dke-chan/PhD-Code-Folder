#
nSim <- as.numeric(commandArgs(trailingOnly = TRUE)[[1]])
nCore <- as.numeric(commandArgs(trailingOnly = TRUE)[[2]])

#
set.seed(2787484)

#
data.frame(seed = sample(0:1e9, nSim), arrayID = rep(1:nCore, length.out = nSim)) |>
  write.csv(file = "./arrayRef.csv", row.names = FALSE)

#
boilerPlate <- data.frame(seed = numeric(0), theta = numeric(0), se = numeric(0), model = character(0), design = character(0))

#
sim.Specs <- LETTERS[1:6]
trap.Specs <- c("3x1_Single", "3x1_Multi")

master.List <- expand.grid(trap.Specs, sim.Specs)

for (i in 1:nrow(master.List)) {
  if (!file.exists(paste0("./results/", master.List[i, 1], "_", master.List[i, 2], ".csv"))) {
    write.csv(boilerPlate, file = paste0("./results/", master.List[i, 1], "_", master.List[i, 2], ".csv"), row.names = FALSE)
  }
}
