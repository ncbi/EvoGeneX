library(ouch)
library(ape)

source("../R/utils.R")

BaseModel = setRefClass("BaseModel",
  fields=list(tree="ouchtree"),
  methods=list(
    plot = function(x) {
      plot(tree)
    },
    setTree = function(newickfile) {
      apetree <- read.tree(newickfile)
      ot <- ape2ouch(apetree, scale=F)
      otd <- as(ot, 'data.frame')
      otd$labels <- as.character(otd$labels)
      otd$labels <- ifelse(otd$labels == "", NA, otd$labels)
      tree <<- with(otd,ouchtree(nodes=nodes,ancestors=ancestors,times=times,labels=labels))
    },
    setRegimes = function(regimefile) {
      regimesdf = read.csv(regimefile, header=T, stringsAsFactors=F)
      nodelabels = tree@nodelabels
      lineages = tree@lineages
      regimesdf$nodeid = sapply(regimesdf$node, function(x) {
        x = unlist(strsplit(x, split="[.]"))
        if (length(x) == 1) {
          return(which(nodelabels == x))
        } else {
          a = which(nodelabels == x[1])
          b = which(nodelabels == x[2])
          lineage1 = rev(lineages[[a]])
          lineage2 = rev(lineages[[b]])
          i = 1
          while (lineage1[i] == lineage2[i]) i = i + 1
          return(lineage1[i-1])
        }
      })
      regimesdf = regimesdf[order(regimesdf$nodeid),]
      regimes <<- data.frame(regimes=factor(regimesdf$regime))
    }
  )
)

