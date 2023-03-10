library(Rcpp); library(RcppArmadillo)

Rcpp.package.skeleton("moveLLR",
                      cpp_files = c("code/gibbons/cpp-Files/ascrStack.cpp", "code/gibbons/cpp-Files/ascrDisperse.cpp"),
                      example_code = FALSE,
                      path = "code/gibbons")
RcppArmadillo.package.skeleton(path = "code/gibbons")

system("mv code/gibbons/anRpackage/src/Makevars code/gibbons/moveLLR/src")
system("mv code/gibbons/anRpackage/src/Makevars.win code/gibbons/moveLLR/src")
system("rm code/gibbons/anRpackage -r")

## Need to update to use reg ex to identify the right lines
DESCRIPTION_lines <- readLines("code/gibbons/moveLLR/DESCRIPTION")
DESCRIPTION_lines[grep("Imports: ", DESCRIPTION_lines)] <- "Imports: Rcpp (>= 0.11.0)"
DESCRIPTION_lines[grep("LinkingTo: ", DESCRIPTION_lines)] <- "LinkingTo: Rcpp, RcppArmadillo"

writeLines(DESCRIPTION_lines, "code/gibbons/moveLLR/DESCRIPTION")
compileAttributes("code/gibbons/moveLLR")

system("R CMD build code/gibbons/moveLLR")
system("mv moveLLR_1.0.tar.gz code/gibbons")

if ("moveLLR" %in% rownames(installed.packages())) {
  try(detach("package:moveLLR", unload = TRUE), silent = TRUE)
  remove.packages("moveLLR")
}

install.packages("code/gibbons/moveLLR_1.0.tar.gz", repos = NULL, type = "source")

system("rm code/gibbons/moveLLR -r")
system("rm code/gibbons/moveLLR_1.0.tar.gz")
