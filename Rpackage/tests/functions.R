# Functions to run Brown and EvoGeneX both with R code and C++ code
run_brown_slow <- function(tree_file, data) {
  brown <- Brown()
  brown$setTree(tree_file)
  results <- brown$fitSlow(data, 0.01, exp_col="exprval")
  return(results)
}

run_brown_fast <- function(tree_file, data) {
  brown <- Brown()
  brown$setTree(tree_file)
  results <- brown$fit(data, 0.01, format="tall", exp_col="exprval")
  return(results)
}

run_evogenex_slow <- function(tree_file, regime_file, data) {
  evog <- EvoGeneX()
  evog$setTree(tree_file)
  evog$setRegimes(regime_file)
  results <- evog$fitSlow(data, 0.01, 0.01)
  return(results)
}

run_evogenex_fast <- function(tree_file, regime_file, data) {
  evog <- EvoGeneX()
  evog$setTree(tree_file)
  evog$setRegimes(regime_file)
  results <- evog$fit(data, format="tall", 0.01, 0.01)
  return(results)
}