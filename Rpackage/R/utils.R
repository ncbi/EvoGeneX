glssoln <- function (a, x, v, tol = sqrt(.Machine$double.eps)) {
  n <- length(x)
  vh <- try(
    chol(v),
    silent=FALSE
  )
  #print("L:")
  #print(vh)
  if (inherits(vh,'try-error')) {
    warning(
      "glssoln: Choleski decomposition of variance-covariance matrix fails",
      call.=FALSE
    )
    y <- rep(NA,ncol(a))
    e <- rep(NA,n)
    dim(y) <- ncol(a)
    dim(e) <- n
  } else {
  #print("X:")
  #print(forwardsolve(vh,a,upper.tri=TRUE,transpose=TRUE))
  #print("y:")
  #print(forwardsolve(vh,x,upper.tri=TRUE,transpose=TRUE))
    s <- svd(forwardsolve(vh,a,upper.tri=TRUE,transpose=TRUE))
    inds <- s$d > tol*max(s$d)
    svals <- s$d[inds,drop=FALSE]
    r <- length(svals)
    svals <-  diag(1/svals,nrow=r,ncol=r)
    y <- (s$v[,inds,drop=FALSE]%*%
            (svals %*%
               t(s$u[,inds,drop=FALSE])))%*%
      forwardsolve(vh,x,upper.tri=TRUE,transpose=TRUE)
    e <- a%*%y-x
    dim(y) <- dim(y)[1]
    dim(e) <- n
  }
  list(
    coeff=y,
    residuals=e
  )
}

sets.of.regimes = function (object, myregimes) {
  lapply(myregimes,function(x)sort(unique(x)))
}

getBeta = function (object, myregimes) {
  nterm <- object@nterm
  nchar <- length(myregimes)
  reg <- sets.of.regimes(object,myregimes)
  nreg <- sapply(reg,length)
  beta <- vector(mode='list',length=nterm)
  #print("lineages")
  #print(object@lineages)
  for (i in seq_len(nterm)) {
    p <- object@lineages[[object@term[i]]]
    np <- length(p)
    beta[[i]] <- vector(mode='list',length=nchar)
    for (n in seq_len(nchar)) {
      beta[[i]][[n]] <- matrix(data=NA,nrow=np,ncol=nreg[n], dimnames = list(NULL, reg[n]$regimes))
      #print(paste("In getBeta for term", i, "char", n, "myreg"))
      #print(myregimes[[n]])
      #print(paste("In getBeta for term", i, "char", n, "reg"))
      #print(reg[[n]])
      #print("p")
      #print(p)
      for (ell in seq_len(nreg[n])) {
        beta[[i]][[n]][,ell] <- ifelse(myregimes[[n]][p]==reg[[n]][ell],1,0)
      }
    }
  }
  beta
}

getBetaCompact = function (object, myregimes) {
  nterm <- object@nterm
  reg <- sets.of.regimes(object,myregimes)
  nreg <- sapply(reg,length)
  beta <- vector(mode='list',length=nterm)
  for (i in seq_len(nterm)) {
    p <- object@lineages[[object@term[i]]]
    beta[[i]] <- as.numeric(myregimes[[1]][p])-1
  }
  as.integer(unlist(beta))
}
