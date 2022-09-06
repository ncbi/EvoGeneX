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


# # Uninstall whatever the current version of EvoGeneX is
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
new_brown_slow <- run_brown_slow(tree_file, data)
new_brown_fast <- run_brown_fast(tree_file, data)
new_evog_slow <- run_evogenex_slow(tree_file, regime_file, data)
new_evog_fast <- run_evogenex_fast(tree_file, regime_file, data)

# Convert results to data frames
nbs <- data.frame("conv"=new_brown_slow$optim.diagn$convergence,
                  "theta"=new_brown_slow$theta,
                  "sigma.sq"=new_brown_slow$sigma.sq,
                  "gamma"=new_brown_slow$gamma,
                  "loglik"=new_brown_slow$loglik)
nbf <- data.frame("conv"=new_brown_fast$optim.diagn$convergence,
                  "theta"=new_brown_fast$theta,
                  "sigma.sq"=new_brown_fast$sigma.sq,
                  "gamma"=new_brown_fast$gamma,
                  "loglik"=new_brown_fast$loglik)
nes <- data.frame("conv"=new_evog_slow$optim.diagn$convergence,
                  "theta"=new_evog_slow$theta,
                  "sigma.sq"=new_evog_slow$sigma.sq,
                  "gamma"=new_evog_slow$gamma,
                  "loglik"=new_evog_slow$loglik)
nef <- data.frame("conv"=new_evog_fast$optim.diagn$convergence,
                  "theta"=new_evog_fast$theta,
                  "sigma.sq"=new_evog_fast$sigma.sq,
                  "gamma"=new_evog_fast$gamma,
                  "loglik"=new_evog_fast$loglik)

# Save results
write.table(nbs, "tests/results/local_brown_slow.csv")
write.table(nbf, "tests/results/local_brown_fast.csv")
write.table(nes, "tests/results/local_evog_slow.csv")
write.table(nef, "tests/results/local_evog_fast.csv")