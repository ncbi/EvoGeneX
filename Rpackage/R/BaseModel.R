#' Base refclass to implement different models
#'
#' @import methods
#' @import ape
#' @import ouch
#' @export BaseModel
#' @exportClass BaseModel

BaseModel = setRefClass("BaseModel",
  fields=list(tree="ouchtree", nbranch="numeric", packed_epochs="numeric"),
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
      nbranch <<- sapply(tree@epochs, length)
      #print(nbranch)
      packed_epochs <<- unlist(tree@epochs)
      #print(packed_epochs)
    }
  )
)

