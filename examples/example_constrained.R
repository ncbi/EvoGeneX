library(EvoGeneX)
library(tidyverse)

wide <- read.csv("../data/HD_M_FBgn0000008.csv", stringsAsFactors = FALSE)
print(wide)

evog <- EvoGeneX()
evog$setTree("../data/drosophila9.newick")
evog$setRegimes("../data/regime_global.csv")

ou_res <- evog$fit(wide, alpha = 0.1, gamma = 0.01)
print(ou_res)

brown <- Brown()
brown$setTree("../data/drosophila9.newick")

brown_res <- brown$fit(wide, gamma = 0.01)
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

print(pvalue)
print(qvalue)

fdr_cutoff <- 0.05

constrained_vs_neutral <- ifelse(qvalue < fdr_cutoff, "constrained", "neutral")

print(constrained_vs_neutral)
