library(EvoGeneX)
library(tidyverse)

wide <- read.csv("../inst/extdata/HD_M_FBgn0000008.csv", stringsAsFactors = FALSE)
cat("\nInput data:\n")
cat("===========\n")
print(wide)

evog <- EvoGeneX()
evog$setTree("../inst/extdata/drosophila9.newick")
evog$setRegimes("../inst/extdata/regime_global.csv")

ou_res <- evog$fit(wide, alpha = 0.1, gamma = 0.01)
cat("\nFit by constrained model:\n")
cat("=========================\n")
print(ou_res)

brown <- Brown()
brown$setTree("../inst/extdata/drosophila9.newick")

brown_res <- brown$fit(wide, gamma = 0.01)
cat("\nFit by neutral model:\n")
cat("=========================\n")
print(brown_res)

# degrees of freedom under different models
ou_dof <- (
  1   # alpha
  + 1 # sigma.sq
  + 1 # theta
  + 1 # gamma
)

brown_dof <- (
  1   # sigma.sq
  + 1 # theta
  + 1 # gamma
)

# loglikelihood ratio test EvoGeneX VS replicated Brown
pvalue <- 1 - pchisq((ou_res$loglik - brown_res$loglik) * 2,
                     (ou_dof - brown_dof))
qvalue <- p.adjust(pvalue, method = "fdr")

cat("\nLoglikelihood ratio test - constrained vs neutral:\n")
cat(paste("P-value:", pvalue, "\n"))
cat(paste("Q-value:", qvalue, "\n"))

fdr_cutoff <- 0.05

constrained_vs_neutral <- ifelse(qvalue < fdr_cutoff, "constrained", "neutral")

cat("\nConsidering all tests, whether constrained or not neutral:\n")
print(constrained_vs_neutral)
