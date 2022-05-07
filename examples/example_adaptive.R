library(EvoGeneX)
library(tidyverse)

wide <- read.csv("../data/HD_M_FBgn0000008.csv", stringsAsFactors = FALSE)
print(wide)

evog <- EvoGeneX()
evog$setTree("../data/drosophila9.newick")
evog$setRegimes("../data/regime_fruitveg.csv")

evog_res <- evog$fit(wide, alpha = 0.1, gamma = 0.01)
print(evog_res)

brown <- Brown()
brown$setTree("../data/drosophila9.newick")

brown_res <- brown$fit(wide, gamma = 0.01)
print(brown_res)
