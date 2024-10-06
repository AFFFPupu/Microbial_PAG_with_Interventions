# requirements.R

# List of required packages
required_packages <- c(
  "dplyr",
  "tidyr",
  "igraph",
  "pcalg",
  "MixedGraphs",
  "deSolve",
  "JuliaCall"
)

# Install any missing packages
installed_packages <- installed.packages()

for (pkg in required_packages) {
  if (!pkg %in% installed_packages[, "Package"]) {
    install.packages(pkg)
  }
}

# If you are using Bioconductor packages (like graph), add this line:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("graph")

# Setup JuliaCall with Julia packages
library(JuliaCall)
julia <- julia_setup()
de <- diffeqr::diffeq_setup()

# Install Julia packages
JuliaCall::julia_library("LinearAlgebra")
JuliaCall::julia_eval('using DifferentialEquations')