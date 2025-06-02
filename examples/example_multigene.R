library(EvoGeneX)
library(tidyverse)

wide <- read.csv("../inst/extdata/HD_M_FBgn0000008.csv", stringsAsFactors = FALSE)
cat("\nInput data (wide format):\n")
cat("=========================\n")
tall <- wide %>% gather("replicate", "exprval", -species)

# pretend we have two different genes
tall <- rbind(tall %>% mutate(symbol = "gene1"),
              tall %>% mutate(symbol = "gene2"))

cat("\nInput data (tall format):\n")
cat("=========================\n")
print(tall)

evog <- EvoGeneX()
evog$setTree("../inst/extdata/drosophila9.newick")
evog$setRegimes("../inst/extdata/regime_global.csv")

brown <- Brown()
brown$setTree("../inst/extdata/drosophila9.newick")

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

fdr_cutoff <- 0.05


process_single_gene <- function(data) {

  ou_res <- evog$fit(data, format = "tall", alpha = 0.1, gamma = 0.01)
  brown_res <- brown$fit(data, format = "tall", gamma = 0.01)

  # loglikelihood ratio test EvoGeneX VS replicated Brownian motion
  pvalue <- 1 - pchisq((ou_res$loglik - brown_res$loglik) * 2,
                       (ou_dof - brown_dof))
}

res <- (
  tall
  %>% group_by(symbol)
  %>% summarize(pvalue = process_single_gene(pick(everything())))
  %>% ungroup()
  %>% mutate(qvalue = p.adjust(pvalue, method = "fdr"))
  %>% mutate( constrained_vs_neutral = ifelse(qvalue < fdr_cutoff,
                                              "constrained", "neutral"))
)
cat("\nResults:\n")
print(res)
