#' EvoGeneX class implementation
#'
#' @import methods
#' @importFrom tidyr gather
#' @export EvoGeneX
#' @exportClass EvoGeneX


EvoGeneX = setRefClass("EvoGeneX",
  contains = "BaseModel",
  fields=list(regimes="data.frame", packed_beta="integer"),
  methods=list(
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
      #print(regimes)
      packed_beta <<- getBetaCompact(tree,regimes)
      #print(packed_beta)
    },
    getWeights = function(nrep, beta, alpha) {
      nt = tree@nterm
      epochs = tree@epochs
      nreg = ncol(beta[[1]][[1]])
      W = matrix(0, nt*nrep, nreg)
      #print("In get weight, initially W:")
      #print(W)
      for (i in 1:nt) {
        ep = epochs[[i]]
        np = length(ep)
        y = vector('numeric', np)
        for (j in 1:np) {
          t = ep[1]-ep[j];
          y[j] = exp(-alpha*t);
        }
        #print("ep")
        #print(ep)
        #print("y")
        #print(y)
        for (j in 1:(np-1)) {
          y[j] = y[j] - y[j+1];
        }
        bp = beta[[i]][[1]]
        #print(paste("BP for terminal", i))
        #print(bp)
        for (r in 1:nreg) {
          for (k in 1:nrep) {
            p = k + (i-1)*nrep
            W[p,r] = sum(y*bp[,r])
          }
        }
      }
      #print("In get weight, at end W:")
      #print(W)
      return(W)
    },
    getCovar = function(nrep, alpha, gamma.sq) {
      nterm = tree@nterm
      bt = tree@branch.times
      v = matrix(0, nterm*nrep, nterm*nrep)
      for (i in 1:nterm) {
        for (k in 1:nrep) {
          for (j in 1:nterm) {
            for (l in 1:nrep) {
              p = k + (i-1)*nrep
              q = l + (j-1)*nrep
              v[p,q] = exp(alpha*(-bt[i,i] - bt[j,j] + 2*bt[i,j]))/(2*alpha)
              if ((i == j) && (k == l)) {
                v[p,q] = v[p,q] + gamma.sq
              }
            }
          }
        }
      }
      return(v)
    },
    computeLogLik = function (nrep, beta, dat, alpha, gamma.sq) {
      n <- length(dat)
      w <- getWeights(nrep, beta, alpha)
      v <- getCovar(nrep, alpha, gamma.sq)

#      print("###################### IN ComputeLogLik ########################")
#      print(paste("n:", n))
#      print(paste("nrep:", nrep))
#      print(paste("alpha:", alpha))
#      print(paste("gamma:", gamma.sq))
#      print("W:")
#      print(w)
#      print("V:")
#      print(v)

      gsol <- try(
        glssoln(w,dat,v),
        silent=FALSE
      )
      if (inherits(gsol,'try-error')) {
        e <- rep(NA,n)
        theta <- rep(NA,ncol(w))
        res <- Inf
      } else {
        e <- gsol$residuals
        theta <- gsol$coeff
        q <- e%*%solve(v,e)
      }
#      print("theta:")
#      print(theta)
#      print("e:")
#      print(e)
#      print("zz:")
#      print(solve(v,e))
#      print("q:")
#      print(q)
      sigma.sq = q[1,1]/n
      det.v = determinant(v, logarithm=T)
      if (det.v$sign != 1) {
        stop("mylogLik error: non-positive determinant",call.=FALSE)
      }
      log.det.v <- det.v$modulus
      res = n*log(2*pi) + n*(1+log(sigma.sq)) + log.det.v
#      print(paste("sigma.sq:", sigma.sq))
#      print(paste("logLik", res))
      list(dev=res, alpha=alpha, gamma.sq=gamma.sq, sigma.sq=sigma.sq, theta=theta)
    },

    fitSlow = function(data, alpha, gamma.sq, lb = 1e-10, ub = 1e+10,...) {
      otd = as(tree, 'data.frame')
      tmp <- merge(otd[c('nodes', 'labels')], data, by.x='labels', by.y='row.names')
      rownames(tmp) <- tmp$nodes
      tmp$nodes <- NULL
      tmp$labels <- NULL
      tmp = tmp[as.character(tree@term),]
      dat = gather(data.frame(t(tmp)))$value
      nrep <- ncol(data)
      beta <- getBeta(tree,regimes)
      optim.diagn <- vector(mode='list',length=0)
      par = c(alpha, gamma.sq)
      opt <- nloptr(par,
                    eval_f = function(par) {
                      computeLogLik(nrep=nrep, beta=beta, dat=dat, alpha=par[1], gamma.sq=par[2])$dev
                    },
                    eval_grad_f=NULL,
                    eval_g_ineq=NULL,
                    eval_g_eq=NULL,
                    eval_jac_g_ineq=NULL,
                    eval_jac_g_eq=NULL,
                    lb=rep(lb, length(par)),
                    ub=rep(ub, length(par)),
                    opts <- list(algorithm="NLOPT_LN_SBPLX",
                                 maxeval=10000,
                                 ftol_rel=.Machine$double.eps^0.5))

      status = opt$status
      message = opt$message

      if (!((status>=1) && (status <= 4))) {
        message("unsuccessful convergence, code ", status, ", see documentation for ", 'nloptr')
        warning("unsuccessful convergence, message ", message)
      }

      optim.diagn <- list(convergence=opt$status,message=opt$message)

      sol <- computeLogLik(nrep=nrep, beta=beta, dat=dat, alpha=opt$solution[1], gamma.sq=opt$solution[2])

      list(optim.diagn=optim.diagn,
        theta=setNames(sol$theta, colnames(beta[[1]][[1]])),
        alpha=sol$alpha,
        sigma.sq=sol$sigma.sq,
        gamma.sq=sol$gamma.sq,
        loglik=-0.5*sol$dev
      )
    },

    fit = function(data, alpha, gamma.sq, lb = 1e-10, ub = 1e+10,...) {
      otd = as(tree, 'data.frame')
      tmp <- merge(otd[c('nodes', 'labels')], data, by.x='labels', by.y='row.names')
      rownames(tmp) <- tmp$nodes
      tmp$nodes <- NULL
      tmp$labels <- NULL
      tmp = tmp[as.character(tree@term),]
      dat = gather(data.frame(t(tmp)))$value
      nrep <- ncol(data)

      opt <- evogenex_fit(dat=dat, nterm=tree@nterm,
                          nreg=nlevels(regimes$regimes),
                          nrep=nrep, nbranch=nbranch, beta=packed_beta,
                          epochs=packed_epochs, bt=tree@branch.times,
                          alpha=alpha, gamma=gamma.sq)

      if (!((opt$status>=1) && (opt$status <= 4))) {
        message("unsuccessful convergence, code ", opt$status, ", see documentation for ", 'nloptr')
        warning("unsuccessful convergence, message ", opt$message)
      }

      optim.diagn <- list(convergence=opt$status,message=opt$message)

      list(optim.diagn=optim.diagn,
        theta=setNames(opt$theta, levels(regimes$regimes)),
        alpha=opt$alpha,
        sigma.sq=opt$sigma.sq,
        gamma.sq=opt$gamma,
        loglik=opt$logLik
      )
    }
  )
)

