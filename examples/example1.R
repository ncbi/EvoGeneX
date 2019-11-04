library(EvoGeneX)

data = read.csv('../pkgdata/HD_M_FBgn0000008.csv', row.names=1, stringsAsFactors=F)
print(data)

evog <- EvoGeneX()
evog$setTree('../pkgdata/drosophila9.newick')
evog$setRegimes('../pkgdata/regime_global.csv')
res = evog$fit(data, alpha=0.1, gamma.sq=0.01)
print(res)

brown <- Brown()
brown$setTree('../pkgdata/drosophila9.newick')
res = brown$fit(data, gamma.sq=0.01)
print(res)
