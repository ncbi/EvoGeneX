library(tidyverse)
library(ape)
library(ouch)
library(reshape2)
library(nloptr)
source("tests/functions.R")

# Installs the local version of EvoGeneX and runs all methods on
# the example drosophila data, saving to tests/results/
# Assumes this is being run from the Rpackage/ directory.

# Test files
file_paths <- get_files()

# Load data
data <- get_data(file_paths["data_file"])

# Uninstall whatever the current version of EvoGeneX is
tryCatch({
    remove.packages("EvoGeneX")
},
error = function(cond) {
    print("EvoGeneX not installed, skipping removal.")
})

# Install the local version
devtools::document()
devtools::build(vignettes=TRUE)
devtools::install(build_vignettes = TRUE)
library(EvoGeneX)

# Run all functions
print("Run with local version.")
run_all(file_paths["tree_file"],
        file_paths["single_regime_file"],
        file_paths["two_regime_file"],
        data,
        "local")
