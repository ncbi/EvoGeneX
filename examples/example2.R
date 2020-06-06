library(EvoGeneX)
library(microbenchmark)

data = read.csv('../data/HD_M_FBgn0000008.csv', row.names=1, stringsAsFactors=F)
print(data)

evog <- EvoGeneX()
evog$setTree('../data/drosophila9.newick')
evog$setRegimes('../data/regime_fruitveg.csv')

#res = evog$fit(data, alpha=0.1, gamma.sq=0.01)
#print(res)
##sink("slow.log", type=c("output"))
#res = evog$fitSlow(data, alpha=0.1, gamma.sq=0.01)
#print(res)
#
#
#mbm <- microbenchmark(
#    "fast" = { res = evog$fit(data, alpha=0.1, gamma.sq=0.01) },
#    "slow" = { res = evog$fitSlow(data, alpha=0.1, gamma.sq=0.01) }
#    )
#print(mbm)
#
#quit()

brown <- Brown()
brown$setTree('../data/drosophila9.newick')
res = brown$fitSlow(data, gamma.sq=0.01)
print(res)
res = brown$fit(data, gamma.sq=0.01)
print(res)
mbm <- microbenchmark(
    "fast" = { res = brown$fit(data, gamma.sq=0.01) },
    "slow" = { res = brown$fitSlow(data, gamma.sq=0.01) }
    )
print(mbm)
