library(EvoGeneX)
library(tidyverse)

wide <- read.csv("../data/HD_M_FBgn0000008.csv", stringsAsFactors = FALSE)
print(wide)

tall <- wide %>% gather("replicate", "exprval", -species)
print(tall)

evog <- EvoGeneX()
evog$setTree("../data/drosophila9.newick")
evog$setRegimes("../data/regime_global.csv")

evog_tall <- evog$fit(tall, format = "tall", alpha = 0.1, gamma = 0.01)
print(evog_tall)

evog_wide <- evog$fit(wide, format = "wide", alpha = 0.1, gamma = 0.01)
print(evog_wide)

brown <- Brown()
brown$setTree("../data/drosophila9.newick")

brown_tall <- brown$fit(wide, format = "wide", gamma = 0.01)
print(brown_tall)

brown_wide <- brown$fit(tall, format = "tall", gamma = 0.01)
print(brown_wide)
