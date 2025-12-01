# Functions to run Brown and EvoGeneX both with R code and C++ code

# Loads the example drosophila data
get_files <- function() {
  return(c(tree_file = "inst/extdata/drosophila9.newick",
           single_regime_file = "inst/extdata/regime_global.csv",
           two_regime_file = "inst/extdata/regime_tworegime.csv",
           data_file = "inst/extdata/HD_M_FBgn0000008.csv"))
}

# Read the data and convert to tall
get_data <- function(data_file) {
  data <- read.csv(data_file)
  data <- data %>% gather("replicate", "exprval", -species)
}

# Run brownian motion R code for a single gene
run_brown_slow <- function(tree_file, data) {
  brown <- Brown()
  brown$setTree(tree_file)
  results <- brown$fitSlow(data, 0.01, exp_col="exprval")
  return(results)
}

# Run brownian motion C++ code for a single gene
run_brown_fast <- function(tree_file, data) {
  brown <- Brown()
  brown$setTree(tree_file)
  results <- brown$fit(data, 0.01, format="tall", exp_col="exprval")
  return(results)
}

# Run EvoGeneX R code for a single gene
run_evogenex_slow <- function(tree_file, regime_file, data) {
  evog <- EvoGeneX()
  evog$setTree(tree_file)
  evog$setRegimes(regime_file)
  results <- evog$fitSlow(data, 0.01, 0.01)
  return(results)
}

# Run EvoGeneX C++ code for a single gene
run_evogenex_fast <- function(tree_file, regime_file, data) {
  evog <- EvoGeneX()
  evog$setTree(tree_file)
  evog$setRegimes(regime_file)
  results <- evog$fit(data, format="tall", 0.01, 0.01)
  return(results)
}

# Check that setting the tree works for both passing a file containing a newick string as well as a newick string directly
run_newick_test <- function(tree_file) {
  evog1 <- EvoGeneX()
  evog1$setTree(tree_file)
  evog2 <- EvoGeneX()
  newick <- readLines(tree_file)
  evog2$setTree(newick)
  t1 <- evog1$tree
  t2 <- evog2$tree
  # Check all fields are equal
  same_tree = TRUE
  if (!identical(t1@nodes, t2@nodes)) {same_tree = FALSE}
  if (same_tree & !identical(t1@ancestors, t2@ancestors)) {same_tree = FALSE}
  if (same_tree & !identical(t1@nodelabels, t2@nodelabels)) {same_tree = FALSE}
  if (same_tree & !identical(t1@times, t2@times)) {same_tree = FALSE}
  if (same_tree & (t1@root != t2@root)) {same_tree = FALSE}
  if (same_tree & (t1@nterm != t2@nterm)) {same_tree = FALSE}
  if (same_tree & !identical(t1@term, t2@term)) {same_tree = FALSE}
  if (same_tree & !identical(t1@anc.numbers, t2@anc.numbers)) {same_tree = FALSE}
  if (same_tree & !identical(t1@lineages, t2@lineages)) {same_tree = FALSE}
  if (same_tree & !identical(t1@epochs, t2@epochs)) {same_tree = FALSE}
  if (same_tree & !identical(t1@branch.times, t2@branch.times)) {same_tree = FALSE}
  if (same_tree & (t1@depth != t2@depth)) {same_tree = FALSE}
  return(same_tree)
}


# Run all methods and write output to tests/results/
run_all <- function(tree_file, single_regime_file, two_regime_file, data, version) {

  print("Running all tests and saving to tests/results/")
  if (!dir.exists(file.path("tests/results/"))){
    print("Directory tests/results/ does not exist. Creating it.")
    dir.create(file.path("tests/results/"))
  }

  print("Run brown slow")
  new_brown_slow <- run_brown_slow(tree_file, data)
  new_brown_slow <- unlist(new_brown_slow)
  write.table(new_brown_slow, paste("tests/results/", version, "_brown_slow.csv", sep=""))

  print("Run brown fast")
  new_brown_fast <- run_brown_fast(tree_file, data)
  new_brown_fast <- unlist(new_brown_fast)
  write.table(new_brown_fast, paste("tests/results/", version, "_brown_fast.csv", sep=""))

  print("Run single regime evogenex slow")
  new_single_evog_slow <- run_evogenex_slow(tree_file, single_regime_file, data)
  new_single_evog_slow <- unlist(new_single_evog_slow)
  write.table(new_single_evog_slow, paste("tests/results/", version, "_single_evog_slow.csv", sep=""))

  print("Run single regime evogenex fast")
  new_single_evog_fast <- run_evogenex_fast(tree_file, single_regime_file, data)
  new_single_evog_fast <- unlist(new_single_evog_fast)
  write.table(new_single_evog_fast, paste("tests/results/", version, "_single_evog_fast.csv", sep=""))

  print("Run two regime evogenex slow")
  new_two_evog_slow <- run_evogenex_slow(tree_file, two_regime_file, data)
  new_two_evog_slow <- unlist(new_two_evog_slow)
  write.table(new_two_evog_slow, paste("tests/results/", version, "_two_evog_slow.csv", sep=""))

  print("Run two regime evogenex fast")
  new_two_evog_fast <- run_evogenex_fast(tree_file, two_regime_file, data)
  new_two_evog_fast <- unlist(new_two_evog_fast)
  write.table(new_two_evog_fast, paste("tests/results/", version, "_two_evog_fast.csv", sep=""))

  print("Run newick file vs string comparison")
  newick_res <- run_newick_test(tree_file)
  write.table(newick_res, paste("tests/results/", version, "_newick.csv", sep=""))
}
