source("../R/EvoGeneX.R")
source("../R/Brown.R")

data = read.csv('../data/HD_M_FBgn0000008.csv', row.names=1, stringsAsFactors=F)
print(data)

evog <- EvoGeneX()
evog$setTree('../data/drosophila9.newick')
evog$setRegimes('../data/regime_fruitveg.csv')
res = evog$fit(data, alpha=0.1, gamma.sq=0.01)
print(res)

brown <- Brown()
brown$setTree('../data/drosophila9.newick')
res = brown$fit(data, gamma.sq=0.01)
print(res)
