library(readxl); library(ascr); library(secr); library(Rcpp); library(RcppArmadillo)

load("code/gibbons/gibbons-prep-data.RData")

x11()
default.par.val <- par()$ask
par(ask = TRUE)
for (i in 1:length(capt.full)) {
  if (nrow(capt.full[[i]]$bincapt) > 0) {
    for (j in 1:nrow(capt.full[[i]]$bincapt)) {
      plot(survey.mask[[i]], type = "n", main = paste0("Session ", i), sub = paste0("Call ", j, ", Animal ", rownames(capt.full[[i]]$bincapt)[j]))
      points(posts.full[[i]], pch = 16, col = "red")

      whichIndices <- capt.full[[i]]$bincapt[j, ] * 1:nrow(posts.full[[i]])
      points(posts.full[[i]][whichIndices > 0, , drop = FALSE], col = "red", cex = 2)

      arrow.length <- 0.03 * min(c(diff(range(survey.mask[[i]][, 1])), diff(range(survey.mask[[i]][, 2]))))

      for (k in whichIndices[whichIndices > 0]) {
        arrows(posts.full[[i]][k, 1], posts.full[[i]][k, 2],
               posts.full[[i]][k, 1] + sin(capt.full[[i]]$bearing[j, k]) * arrow.length,
               posts.full[[i]][k, 2] + cos(capt.full[[i]]$bearing[j, k]) * arrow.length,
               0.1, lwd = 2, col = "red")
      }
    }
  }
}
par(ask = default.par.val)
readline()
dev.off()
