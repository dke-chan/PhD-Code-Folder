library(Rcpp); library(RcppArmadillo)

Rcpp.package.skeleton("moveLLR", cpp_files = c("./ascrStack.cpp", "./ascrDisperse.cpp"), example_code = FALSE)
RcppArmadillo.package.skeleton()

system("mv ./anRpackage/src/Makevars ./moveLLR/src")
system("mv ./anRpackage/src/Makevars.win ./moveLLR/src")
system("rm ./anRpackage -r")

## Need to update to use reg ex to identify the right lines
DESCRIPTION_lines <- readLines("./moveLLR/DESCRIPTION")
DESCRIPTION_lines[grep("Imports: ", DESCRIPTION_lines)] <- "Imports: Rcpp (>= 0.11.0)"
DESCRIPTION_lines[grep("LinkingTo: ", DESCRIPTION_lines)] <- "LinkingTo: Rcpp, RcppArmadillo"

writeLines(DESCRIPTION_lines, "./moveLLR/DESCRIPTION")
compileAttributes("./moveLLR")

system("R CMD build ./moveLLR")

if ("moveLLR" %in% rownames(installed.packages())) {
  try(detach("package:moveLLR", unload = TRUE), silent = TRUE)
  remove.packages("moveLLR")
}

install.packages("./moveLLR_1.0.tar.gz", repos = NULL, type = "source")

system("rm ./moveLLR -r")
system("rm ./moveLLR_1.0.tar.gz")
