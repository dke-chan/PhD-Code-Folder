interactive <- FALSE
source("code/gibbons/helper-Scripts/analysis-konbini.r")

# Figure 4.3 --------------------------------------------------------------
pdf("./figures/04-figure03.pdf", height = 6, width = 6)
# layout(matrix(c(1, 2), 1, 2))

listeningPostsSingle <- expand.grid(x = c(-1000, 0, 1000), y = c(0))
listeningPostsMulti.list <- list(expand.grid(x = c(-11000, -10000, -9000), y = 7500),
                                 expand.grid(x = c(11000, 10000, 9000), y = 7500),
                                 expand.grid(x = c(-1000, 0, 1000), y = c(0)),
                                 expand.grid(x = c(-11000, -10000, -9000), y = -7500),
                                 expand.grid(x = c(11000, 10000, 9000), y = -7500))
listeningPostsMulti <- do.call(what = "rbind", listeningPostsMulti.list)

# par(mai = c(0.55, 0.55, 0.1, 0.1))
#
# plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
#      xlim = range(listeningPostsMulti$x) + c(-1, 1) * 5000 ,
#      ylim = range(listeningPostsMulti$y) + c(-1, 1) * 5000)
#
# axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0)); mtext("X", side = 1, line = 1.55, font = 2)
# axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.5, 0)); mtext("Y", side = 2, line = 1.4, font = 2)
#
# points(listeningPostsSingle, pch = 16, cex = 1.1)
#
# mask.1 <- ascr::create.mask(listeningPostsSingle, buffer = 5000)
#   chull(mask.1) |>
#   {\(x) {lines(mask.1[c(x, x[1]), ], lty = 5, lwd = 2, col = "tomato")}}()
#
# box()

par(mai = c(0.55, 0.55, 0.1, 0.1))


plot(0, type = "n", axes = FALSE, xlab = "", ylab = "",
     xlim = range(listeningPostsMulti$x) + c(-1, 1) * 5000,
     ylim = range(listeningPostsMulti$y) + c(-1, 1) * 5000)

lapply(listeningPostsMulti.list, \(a) {
  mask.a <- ascr::create.mask(a, buffer = 5000)
  chull(mask.a) |>
    {\(x) {lines(mask.a[c(x, x[1]), ], lty = 5, lwd = 2, col = "tomato")}}()
})

axis(1, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.6, 0)); mtext("X", side = 1, line = 1.55, font = 2, cex = 1.1)
axis(2, tck = -0.01, cex.axis = 1.05, mgp = c(1, 0.5, 0)); mtext("Y", side = 2, line = 1.4, font = 2, cex = 1.1)

points(listeningPostsMulti, pch = 16, cex = 1.1)

box()

dev.off()

# Data cleaning -----------------------------------------------------------

# files <- c("3x1_Single_A.csv", "3x1_Single_B.csv", "3x1_Single_C.csv", "3x1_Single_D.csv", "3x1_Single_E.csv", "3x1_Single_F.csv",
           # "3x1_Multi_A.csv", "3x1_Multi_B.csv",  "3x1_Multi_C.csv",  "3x1_Multi_D.csv", "3x1_Multi_E.csv", "3x1_Multi_F.csv")
files <- c("3x1_Multi_A.csv", "3x1_Multi_C.csv",  "3x1_Multi_E.csv")
binomialCallProcess <- rep(TRUE, 12) |>
  setNames(nm = files)

masterList <- data.frame(seed = numeric(0), theta = numeric(0), se = numeric(0), model = character(0), design = character(0), callProcess = character(0), coefName = character(0))

for (index in files) {
  simData_i <- read.csv(paste0("code/gibbons/sim-Studies/results/", index))
  seeds <- unique(simData_i$seed)

  if (index == "3x1_Single_A.csv" | index == "3x1_Multi_E.csv") {
    doubleBlocks <- sapply(seeds, \(x) {
      nrow(simData_i[which(simData_i$seed == x), ]) > 17
    })

    ##
    for (seed_i in seeds[which(doubleBlocks)]) {
      if (interactive) {
        cat("Sim set:", index, "\n")
        subset(simData_i, seed == seed_i) |> print()
        readline()
      }
      simData_i <- simData_i[-c(which(simData_i$seed == seed_i)[1:17]), ]
    }
  }

  if (binomialCallProcess[index]) {
    simData_i$callProcess <- "Binomial"
    simData_i$coefName <- ""
    for (seed_i in seeds) {
      simData_i$coefName[which(simData_i$seed == seed_i)] <- c("D_c", "lambda0", "sigma", "kappa", "ell_D", "ell_lambda0", "ell_sigma", "logit_pb", "ell_kappa", "cov_ell_D_logit_pb", "ell_D", "ell_lambda0", "ell_sigma", "logit_pb", "ell_sigmam", "ell_kappa", "cov_ell_D_logit_pb")

      for (model_i in c("ascr", "ascrStack", "ascrDisperse")) {
        model.indices_i <- which(simData_i$seed == seed_i & simData_i$model == model_i)
        if (all(simData_i$theta[model.indices_i] == 0, na.rm = TRUE)) {
          simData_i$theta[model.indices_i] <- NA
        }

        if (all(simData_i$se[model.indices_i] == 0, na.rm = TRUE)) {
          simData_i$se[model.indices_i] <- NA
        }
      }
    }
  }

  ##
  masterList <- rbind(masterList, simData_i)
}

# Report Setup ------------------------------------------------------------

Design_A_Theta <- c(D_c = 0.005 * 2/3 * 3, nc = 3, D = 0.005, lambda0 = 5, sigma = 750, pb = 2/3, sigmam = sqrt((30 * 10000) / pi) / 2.4477, kappa = 10)
Design_B_Theta <- c(D_c = 0.005 * 2/3 * 5, nc = 5, D = 0.005, lambda0 = 5, sigma = 750, pb = 2/3, sigmam = sqrt((30 * 10000) / pi) / 2.4477, kappa = 10)
Design_C_Theta <- c(D_c = 0.005 * 2/3 * 3, nc = 3, D = 0.005, lambda0 = 5, sigma = 750, pb = 2/3, sigmam = sqrt((30 * 10000) / pi), kappa = 10)
Design_D_Theta <- c(D_c = 0.005 * 2/3 * 5, nc = 5, D = 0.005, lambda0 = 5, sigma = 750, pb = 2/3, sigmam = sqrt((30 * 10000) / pi), kappa = 10)
Design_E_Theta <- c(D_c = 0.005 * 2/3 * 3, nc = 3, D = 0.005, lambda0 = 5, sigma = 750, pb = 2/3, sigmam = 1000, kappa = 10)
Design_F_Theta <- c(D_c = 0.005 * 2/3 * 5, nc = 5, D = 0.005, lambda0 = 5, sigma = 750, pb = 2/3, sigmam = 1000, kappa = 10)

Sim.Theta <- list(
  "3x1_Single_A" = Design_A_Theta, "3x1_Multi_A" = Design_A_Theta, "3x1_Single_B" = Design_B_Theta, "3x1_Multi_B" = Design_B_Theta,
  "3x1_Single_C" = Design_C_Theta, "3x1_Multi_C" = Design_C_Theta, "3x1_Single_D" = Design_D_Theta, "3x1_Multi_D" = Design_D_Theta,
  "3x1_Single_E" = Design_E_Theta, "3x1_Multi_E" = Design_E_Theta, "3x1_Single_F" = Design_F_Theta, "3x1_Multi_F" = Design_F_Theta
)

design.labels <- paste0(
  rep(c("C", "S", "D"), times = length(Sim.Theta)),
  sep = "-",
  rep(c("S", "M"), each = 3))

design.cols <- RColorBrewer::brewer.pal(7, "Pastel1")[2:7]

# -------------------------------------------------------------------------

masterList$model.formal <- with(masterList, ifelse(model == "ascr", "ASCR-CD", ifelse(model == "ascrStack", "ASCR-AD-S", "ASCR-AD-D"))) |>
  factor(levels = c("ASCR-CD", "ASCR-AD-S", "ASCR-AD-D"), labels = c("Standard", "Call Stacking", "Call Dispersal"))
masterList$design.formal <- with(masterList, ifelse(grepl("Single", design), "Single Group", "Multiple Groups")) |>
  factor(levels = c("Single Group", "Multiple Groups"))
masterList$nc.sim <- with(masterList, ifelse(grepl("*_A", design) | grepl("*_C", design) | grepl("*_E", design), "3", "5")) |>
  factor(levels = c("3", "5"))
masterList$sigma.sim <- with(masterList, ifelse(grepl("*_A", design) | grepl("*_B", design), "126",
                                             ifelse(grepl("*_C", design) | grepl("*_D", design), "309",
                                                    "1000"))) |>
  factor(levels = c("126", "309", "1000"))

levels(masterList$nc.sim) <- c(expression("n"["c"]==3), expression("n"["c"]==5))
# levels(masterList$sigma.sim) <- c(expression(sigma["m"]==126), expression(sigma["m"]==309), expression(sigma["m"]==1000))

MyPalette <- palette.colors(n = 3, palette = "Accent")
# names(MyPalette) <- c("ASCR-CD", "ASCR-AD-S", "ASCR-AD-D")
names(MyPalette) <- c("Standard", "Call Stacking", "Call Dispersal")

# -------------------------------------------------------------------------

masterList.ggplot <- masterList

refGrid <- expand.grid(model = unique(masterList$model), design = unique(masterList$design), stringsAsFactors = FALSE)

processed <- split(refGrid, ~ model + design) |>
  lapply(\(x) {
    y <- thetaExtract(masterList, x$model, x$design, Sim.Theta[[x$design]], FALSE, include.derived = TRUE)
    lapply(y, \(z) z$estimates)
  })

masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "ell_D", "D", masterList.ggplot$coefName)
masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "ell_lambda0", "lambda0", masterList.ggplot$coefName)
masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "ell_sigma", "sigma", masterList.ggplot$coefName)
masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "logit_pb", "pb", masterList.ggplot$coefName)
masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "ell_kappa", "kappa", masterList.ggplot$coefName)
masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "ell_sigmam", "sigmam", masterList.ggplot$coefName)
masterList.ggplot$coefName <- ifelse(masterList.ggplot$coefName == "cov_ell_D_logit_pb", "D_c", masterList.ggplot$coefName)

for (i in names(processed)) {
  model.name = gsub("(ascr|ascrDisperse|ascrStack)\\.[0-9]x[0-9]_(Multi|Single)_[A-Z]", "\\1", i)
  design.name = gsub("(ascr|ascrDisperse|ascrStack)\\.([0-9]x[0-9]_(Multi|Single)_[A-Z])", "\\2", i)

  for (j in names(processed[[i]])) {

    bool.check = which(masterList.ggplot$design == design.name & masterList.ggplot$model == model.name & masterList.ggplot$coefName == j)
    if (length(bool.check) > 0) {
      masterList.ggplot$theta[bool.check] <- processed[[i]][[j]]
      print(j)
    }
  }

  cat(model.name, "\n")
  cat(design.name, "\n")
}

# -------------------------------------------------------------------------

library(ggplot2)

# ggplot.text = expand.grid(nc.sim = levels(masterList$nc.sim),
            # sigma.sim = levels(masterList$sigma.sim)) |>
  # cbind(data.frame(x = -Inf, y = Inf, model.formal = "ASCR-AD-D"), label = "𝜽", labelno = paste0("phantom(0)[", 1:6, "]"))

##
ggplot(data = subset(masterList.ggplot, coefName == "D" ),
       aes(x = sigma.sim, y = theta * 100, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  # facet_grid(sigma.sim ~ nc.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette[2:3]) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat(D)*" "*bgroup("(","animals km"^-2,")")), fill = "Model") +
  geom_hline(yintercept = 0.5, lty = 5, col = "#f0027f") +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN1.pdf", width = 8, height = 5.5)
# 2 * 1000 * 12 - subset(masterList.ggplot, coefName == "D" & theta <= 0.03) |> nrow()
subset(masterList.ggplot, coefName == "D")$theta |> is.na() |> sum()

ggplot(data = subset(masterList.ggplot, coefName == "lambda0" & theta < 50),
       aes(x = sigma.sim, y = theta, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  # facet_grid(sigma.sim ~ nc.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat(lambda)[0]), fill = "Model") +
  geom_hline(yintercept = 5, lty = 5, col = "#f0027f") +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN2.pdf", width = 8, height = 5.5)
subset(masterList.ggplot, coefName == "lambda0" & theta >= 50) |> nrow()
subset(masterList.ggplot, coefName == "lambda0")$theta |> is.na() |> sum()

##
ggplot(data = subset(masterList.ggplot, coefName == "sigma"),
       aes(x = sigma.sim, y = theta, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  # facet_grid(sigma.sim ~ nc.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat(sigma)*" (m)"), fill = "Model") +
  geom_hline(yintercept = 750, lty = 5, col = "#f0027f") +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN3.pdf", width = 8, height = 5.5)
# 2 * 1000 * 18 - subset(masterList.ggplot, coefName == "sigma") |> nrow()
subset(masterList.ggplot, coefName == "sigma")$theta |> is.na() |> sum()

##
ggplot(data = subset(masterList.ggplot, coefName == "kappa"),
       aes(x = sigma.sim, y = theta, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  # facet_grid(sigma.sim ~ nc.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat(kappa)), fill = "Model") +
  geom_hline(yintercept = 10, lty = 5, col = "#f0027f") +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN4.pdf", width = 8, height = 5.5)
# 2 * 1000 * 18 - subset(masterList.ggplot, coefName == "kappa" & theta < 100) |> nrow()
subset(masterList.ggplot, coefName == "kappa")$theta |> is.na() |> sum()

##
ggplot(data = subset(masterList.ggplot, coefName == "pb"),
       aes(x = sigma.sim, y = theta, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  # facet_grid(sigma.sim ~ nc.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette[2:3]) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat("p")["b"]), fill = "Model") +
  geom_hline(yintercept = 2/3, lty = 5, col = "#f0027f") +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN5.pdf", width = 8, height = 5.5)
subset(masterList.ggplot, coefName == "pb")$theta |> is.na() |> sum()

##
hline.dummy = expand.grid(nc.sim = levels(masterList$nc.sim), sigma.sim = levels(masterList$sigma.sim))
hline.dummy$theta = c(sqrt((30 * 10000) / pi) / 2.4477, sqrt((30 * 10000) / pi) / 2.4477, sqrt((30 * 10000) / pi), sqrt((30 * 10000) / pi), 1000, 1000)

ggplot(data = subset(masterList.ggplot, coefName == "sigmam"),
       aes(x = sigma.sim, y = theta, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  geom_point(aes(x = c(1,1, 2,2, 3, 3), y = theta), data = hline.dummy, size = 3, pch = 23, fill = "#f0027f", inherit.aes = FALSE) +
  # facet_grid( ~ sigma.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette[3]) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat(sigma)["m"]*" (m)"), fill = "Model", pch = 21) +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN6.pdf", width = 8, height = 5.5)
# 2 * 1000 * 3 - subset(masterList.ggplot, coefName == "sigmam" & theta < 2500) |> nrow()
subset(masterList.ggplot, coefName == "sigmam")$theta |> is.na() |> sum()

##
# hline2.dummy = expand.grid(nc.sim = levels(masterList$nc.sim), sigma.sim = levels(masterList$sigma.sim))
# hline2.dummy$theta = c(0.5 * 3 * 2/3, 0.5 * 5 * 2/3, 0.5 * 3 * 2/3, 0.5 * 5 * 2/3, 0.5 * 3 * 2/3, 0.5 * 5 * 2/3)

ggplot(data = subset(masterList.ggplot, coefName == "D_c"),
       aes(x = sigma.sim, y = theta * 100 / 3, fill = model.formal)) +
  geom_boxplot(outlier.alpha = 0.3, width = 0.5, position = position_dodge(0.75)) +
  # geom_text(aes(x = x, y = y, label = label), data = ggplot.text, hjust = -0.75, vjust = 1.2, size = 6) +
  # geom_text(aes(x = x, y = y, label = labelno), data = ggplot.text, hjust = -0.75, vjust = 1.7, size = 4, parse = TRUE) +
  geom_hline(yintercept = 0.5 * 2/3, lty = 5, col = "#f0027f") +
  # facet_grid(sigma.sim ~ nc.sim, labeller = label_parsed) +
  theme_bw() +
  scale_fill_manual(values = MyPalette) +
  labs(x = expression(sigma["m"]*" (m)"), y = expression(hat("D")["c"]*" "*bgroup("(","calls km"^-2*" day"^-1,")")), fill = "Model") +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.background = element_rect(fill = "lemonchiffon"),
        text = element_text(size = 14), legend.margin = margin(-10, 0, 0, 0))
ggsave("figures/04-figureN7.pdf", width = 8, height = 5.5)
# 2 * 1000 * 18 - subset(masterList.ggplot, coefName == "D_c" & theta < 0.1) |> nrow()
subset(masterList.ggplot, coefName == "D_c")$theta |> is.na() |> sum()

# Table Prep --------------------------------------------------------------

toPrint <- thetaResults(masterList, whichDesigns = names(Sim.Theta)[c(2, 6, 10)], truths = Sim.Theta, exclude.outliers = TRUE,
                        derived = TRUE, clean = FALSE)

toPrint$prettyDesign <- paste0(
  gsub("3x1_(Single|Multi)_([A-Z])", "\\2", toPrint$Design),
  sep = " -- ",
  gsub("3x1_(Single|Multi)_[A-Z]", "\\1", toPrint$Design)
)

toPrint$prettyModel <- ""
toPrint$prettyModel[toPrint$Model == "ascr"] <- "ASCR-CD"
toPrint$prettyModel[toPrint$Model == "ascrStack"] <- "ASCR-AD-S"
toPrint$prettyModel[toPrint$Model == "ascrDisperse"] <- "ASCR-AD-D"

toPrint$prettyTheta <- ""
toPrint$prettyTheta[toPrint$Theta == "D_c"] <- "$D_c$"
toPrint$prettyTheta[toPrint$Theta == "lambda0"] <- "$\\lambda_0$"
toPrint$prettyTheta[toPrint$Theta == "sigma"] <- "$\\sigma$"
toPrint$prettyTheta[toPrint$Theta == "kappa"] <- "$\\kappa$"
toPrint$prettyTheta[toPrint$Theta == "D"] <- "$D$"
toPrint$prettyTheta[toPrint$Theta == "pb"] <- "$p_b$"
toPrint$prettyTheta[toPrint$Theta == "sigmam"] <- "$\\sigma_m$"

library(xtable)

tablePrintr <- function (param, muSprintf = "%.2f") {
  x <- subset(toPrint, Theta == param)[, c("prettyDesign", "prettyModel", "prettyTheta",  "Mean", "pctBias", "CV", "waldCapture")]

  for (i in unique(x$prettyDesign)) {
    y <- which(x$prettyDesign == i, arr.ind	= TRUE)
    x$prettyDesign[y[-1]] <- ""
  }

  x.size <- length(y)

  if (param == "D" | param == "D_c") {
    x$Mean <- x$Mean * 100
  }

  x$Mean <- sprintf(muSprintf, x$Mean)
  x$pctBias <- sprintf("%.1f", x$pctBias)
  x$CV <- sprintf("%.1f", x$CV)
  x$waldCapture <- sprintf("%.1f", x$waldCapture)

  colnames(x) <- c("\\textbf{Scenario}", "\\textbf{Model}", "\\textbf{Parameter}",  "\\textbf{Mean}", "\\textbf{Bias (\\%)}", "\\textbf{CV (\\%)}", "\\textbf{Capture (\\%)}")

  xtable(x, "I'm a caption", "I'm a label") |>
    print(include.rownames = FALSE,
          hline.after = c(-1, seq(0, nrow(x), by = x.size)),
          sanitize.text.function = \(a) {a},
          sanitize.colnames.function = \(a) {a},
          booktabs = TRUE)
}

# Tables ------------------------------------------------------------------

tablePrintr("D")
tablePrintr("D_c")

tablePrintr("lambda0", "%.1f")
tablePrintr("sigma", "%.1f")
tablePrintr("kappa", "%.1f")
tablePrintr("pb")
tablePrintr("sigmam", "%.1f")
