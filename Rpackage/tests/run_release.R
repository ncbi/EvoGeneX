library(tidyverse)
library(ape)
library(ouch)
library(reshape2)
library(nloptr)

source("tests/functions.R")

# Test files
tree_file <- "inst/extdata/drosophila9.newick"
regime_file <- "inst/extdata/regime_global.csv"
data_file <- "inst/extdata/HD_M_FBgn0000008.csv"

# Load data
data <- read.csv(data_file)
data <- data %>% gather("replicate", "exprval", -species)

# Uninstall whatever the current version of EvoGeneX is
tryCatch({
    remove.packages("EvoGeneX")
},
error = function(cond) {
    print("EvoGeneX not installed, skipping removal.")
})

# Install current EvoGeneX version
print("Install current EvoGeneX version.")
devtools::install_github("ncbi/EvoGeneX", subdir="Rpackage", build_vignettes=TRUE)
library(EvoGeneX)

# Run all functions
print("Running with current EvoGeneX version.")
curr_brown_slow <- run_brown_slow(tree_file, data)
curr_brown_fast <- run_brown_fast(tree_file, data)
curr_evog_slow <- run_evogenex_slow(tree_file, regime_file, data)
curr_evog_fast <- run_evogenex_fast(tree_file, regime_file, data)

cbs <- data.frame("conv"=curr_brown_slow$optim.diagn$convergence,
                  "theta"=curr_brown_slow$theta,
                  "sigma.sq"=curr_brown_slow$sigma.sq,
                  "gamma"=curr_brown_slow$gamma,
                  "loglik"=curr_brown_slow$loglik)
cbf <- data.frame("conv"=curr_brown_fast$optim.diagn$convergence,
                  "theta"=curr_brown_fast$theta,
                  "sigma.sq"=curr_brown_fast$sigma.sq,
                  "gamma"=curr_brown_fast$gamma,
                  "loglik"=curr_brown_fast$loglik)
ces <- data.frame("conv"=curr_evog_slow$optim.diagn$convergence,
                  "theta"=curr_evog_slow$theta,
                  "sigma.sq"=curr_evog_slow$sigma.sq,
                  "gamma"=curr_evog_slow$gamma,
                  "loglik"=curr_evog_slow$loglik)
cef <- data.frame("conv"=curr_evog_fast$optim.diagn$convergence,
                  "theta"=curr_evog_fast$theta,
                  "sigma.sq"=curr_evog_fast$sigma.sq,
                  "gamma"=curr_evog_fast$gamma,
                  "loglik"=curr_evog_fast$loglik)

write.table(cbs, "tests/results/curr_brown_slow.csv")
write.table(cbf, "tests/results/curr_brown_fast.csv")
write.table(ces, "tests/results/curr_evog_slow.csv")
write.table(cef, "tests/results/curr_evog_fast.csv")