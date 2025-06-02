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

evog2 <- EvoGeneX()
evog2$setTree("../inst/extdata/drosophila9.newick")
evog2$setRegimes("../inst/extdata/regime_tworegime.csv")

ou2_res <- evog2$fit(wide, alpha = 0.1, gamma = 0.01)
cat("\nFit by adaptive model:\n")
cat("=========================\n")
print(ou2_res)

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

ou2_dof <- (
  1   # alpha
  + 1 # sigma.sq
  + 2 # theta
  + 1 # gamma
)

brown_dof <- (
  1   # sigma.sq
  + 1 # theta
  + 1 # gamma
)

# loglikelihood ratio test OU2 VS replicated Brownian Motion
ou2_vs_bm_pvalue <- 1 - pchisq((ou2_res$loglik - brown_res$loglik) * 2,
                               (ou2_dof - brown_dof))
ou2_vs_bm_qvalue <- p.adjust(ou2_vs_bm_pvalue, method = "fdr")
cat("\nQ-value of loglikelihood ratio test - adaptive vs neutral:\n")
print(ou2_vs_bm_qvalue)

# loglikelihood ratio test OU2 VS OU
ou2_vs_ou_pvalue <- 1 - pchisq((ou2_res$loglik - ou_res$loglik) * 2,
                               (ou2_dof - ou_dof))
ou2_vs_ou_qvalue <- p.adjust(ou2_vs_ou_pvalue, method = "fdr")

cat("\nQ-value of loglikelihood ratio test - adaptive vs constrained:\n")
print(ou2_vs_ou_qvalue)

fdr_cutoff <- 0.05

# gene is adaptive if both OU2vsBM and OU2vsOU tests are successfull
is_adaptive <- ifelse(max(ou2_vs_bm_qvalue, ou2_vs_ou_qvalue) < fdr_cutoff,
                      "adaptive", "not-adaptive")

cat("\nConsidering all tests, whether adaptive or not adaptive:\n")
print(is_adaptive)
