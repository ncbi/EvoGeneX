#' EvoGeneX class implementation
#'
#' @import methods
#' @import reshape2
#' @import tidyr
#' @export EvoGeneX
#' @exportClass EvoGeneX


EvoGeneX <- setRefClass("EvoGeneX",
  contains = "BaseModel",
  fields = list(regimes = "data.frame", packed_beta = "integer"),
  methods = list(
    setRegimes = function(regimefile) {
      regimesdf <- read.csv(regimefile, stringsAsFactors = FALSE)
      nodelabels <- tree@nodelabels
      lineages <- tree@lineages
      get_mra <- function(lab1, lab2) {
        if (lab2 %in% c("", NA)) {
          # use node label directly
          return(which(nodelabels == lab1))
        } else {
          # find most recent ancestor of two leaves given by the labels
          id1 <- which(nodelabels == lab1)
          id2 <- which(nodelabels == lab2)
          lineage1 <- rev(lineages[[id1]])
          lineage2 <- rev(lineages[[id2]])
          i <- 1
          while (lineage1[i] == lineage2[i]) i <- i + 1
          return(lineage1[[i-1]])
        }
      }
      regimesdf <- rowwise(regimesdf) %>% mutate(nodeid = get_mra(node, node2))
      if (length(setdiff(tree@nodes, as.character(regimesdf$nodeid))) > 0) {
        stop("There are tree nodes (",
          setdiff(tree@nodes, as.character(regimesdf$nodeid)),
          ") for which regime has not been set!!")
      }
      regimesdf <- regimesdf[order(regimesdf$nodeid),]
      regimes <<- data.frame(regimes=factor(regimesdf$regime))
      packed_beta <<- getBetaCompact(tree,regimes)
    },
    getWeights = function(nrep, beta, alpha) {
      nt = tree@nterm
      epochs = tree@epochs
      nreg = ncol(beta[[1]][[1]])
      W = matrix(0, sum(nrep$n), nreg)
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
          for (k in 1:nrep[i,]$n) {
            p = k + sum(nrep[1:i-1,]$n)
            W[p,r] = sum(y*bp[,r])
          }
        }
      }
      #print("In get weight, at end W:")
      #print(W)
      return(W)
    },
    getCovar = function(nrep, alpha, gamma) {
      nterm = tree@nterm
      bt = tree@branch.times
      sum_nrep = sum(nrep$n)
      v = matrix(0, sum_nrep, sum_nrep)
      for (i in 1:nterm) {
        for (k in 1:nrep[i,]$n) {
          for (j in 1:nterm) {
            for (l in 1:nrep[j,]$n) {
              p = k + sum(nrep[1:i-1,]$n)
              q = l + sum(nrep[1:j-1,]$n)
              part1 = exp(alpha*(-bt[i,i] - bt[j,j] + 2*bt[i,j]));
              part2 = (1-exp(-2*alpha*bt[i,j]));
              v[p,q] = part1*part2/(2*alpha);
              if ((i == j) && (k == l)) {
                v[p,q] = v[p,q] + gamma
              }
            }
          }
        }
      }
      return(v)
    },
    computeLogLik = function (nrep, beta, dat, alpha, gamma) {
      n <- length(dat)
      w <- getWeights(nrep, beta, alpha)
      v <- getCovar(nrep, alpha, gamma)

#      print("###################### IN ComputeLogLik ########################")
#      print(paste("n:", n))
#      print(paste("nrep:", nrep))
#      print(paste("alpha:", alpha))
#      print(paste("gamma:", gamma))
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
        stop("myloglik error: non-positive determinant",call.=FALSE)
      }
      log.det.v <- det.v$modulus
      res = n*log(2*pi) + n*(1+log(sigma.sq)) + log.det.v
#      print(paste("sigma.sq:", sigma.sq))
#      print(paste("loglik", res))
      list(dev=res, alpha=alpha, gamma=gamma, sigma.sq=sigma.sq, theta=theta)
    },



    fitMCMC = function(data, alpha, gamma, lb = 1e-10, ub = 1e+10,...) {
      otd = as(tree, 'data.frame')
      tmp <- merge(otd[c('nodes', 'labels')], data, by.x='labels', by.y='row.names')
      rownames(tmp) <- tmp$nodes
      tmp$nodes <- NULL
      tmp$labels <- NULL
      tmp = tmp[as.character(tree@term),]
      dat = tidyr::gather(data.frame(t(tmp)))$value
      nrep <- ncol(data)
      beta <- getBeta(tree,regimes)
      optim.diagn <- vector(mode='list',length=0)
      par = c(alpha, gamma)
      opt <- nloptr(par,
                    eval_f = function(par) {
                      computeLogLik(nrep=nrep, beta=beta, dat=dat, alpha=par[1], gamma=par[2])$dev
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

      sol <- computeLogLik(nrep=nrep, beta=beta, dat=dat, alpha=opt$solution[1], gamma=opt$solution[2])

      list(optim.diagn=optim.diagn,
        theta=setNames(sol$theta, colnames(beta[[1]][[1]])),
        alpha=sol$alpha,
        sigma.sq=sol$sigma.sq,
        gamma=sol$gamma,
        loglik=-0.5*sol$dev
      )
    },


    fitSlow = function(data, alpha, gamma,
                   species_col = 'species', replicate_col = 'replicate', exprval_col = 'exprval',
                   lb = 1e-10, ub = 1e+10,...) {
      data = data[complete.cases(data),]
      dat = data[c(species_col, replicate_col, exprval_col)]
      names(dat) = c('species', 'replicate', 'exprval')
      replicates = unique(dat$replicate)
      nrep=dat %>% group_by(species) %>% tally()
      dat = dcast(dat, species~replicate, value.var='exprval')
      otd = as(tree, 'data.frame')
      tmp <- merge(otd, data.frame(dat), by.x='labels', by.y='species', all=TRUE)
      # merging destroys index of dataframe
      rownames(tmp) <- tmp$nodes
      tmp = tmp[as.character(tree@term), replicates, drop=FALSE]
      dat = gather(data.frame(t(tmp)))$value
      dat = dat[!is.na(dat)]
      beta <- getBeta(tree,regimes)
      optim.diagn <- vector(mode='list',length=0)
      par = c(alpha, gamma)
      opt <- nloptr(par,
                    eval_f = function(par) {
                      computeLogLik(nrep=nrep, beta=beta, dat=dat, alpha=par[1], gamma=par[2])$dev
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

      sol <- computeLogLik(nrep=nrep, beta=beta, dat=dat, alpha=opt$solution[1], gamma=opt$solution[2])

      list(optim.diagn=optim.diagn,
        theta=setNames(sol$theta, colnames(beta[[1]][[1]])),
        alpha=sol$alpha,
        sigma.sq=sol$sigma.sq,
        gamma=sol$gamma,
        loglik=-0.5*sol$dev
      )
    },

    # input data in melted form with three columns species, replicate, exp
    # where each row of the form sp,rep,ex represent the expression value
    # ex of species sp and replicate rep, example:
    # dmel,R1,123.456 and so on

    fit = function(data, alpha, gamma,
                   format = "wide",
                   species_col = 'species',
                   replicate_col = 'replicate',
                   exprval_col = 'exprval',
                   lb = 1e-10,
                   ub = 1e+10,
                   ...) {
      dat_nrep <- prepare_replicated_data(data,
                                     format,
                                     species_col,
                                     replicate_col,
                                     exprval_col,
                                     tree)

      dat <- dat_nrep$dat
      nrep <- dat_nrep$nrep

      opt <- evogenex_fit(dat = dat,
                          nterm = tree@nterm,
                          nreg = nlevels(regimes$regimes),
                          nrep = nrep,
                          nbranch = nbranch,
                          beta = packed_beta,
                          epochs = packed_epochs,
                          bt = tree@branch.times,
                          alpha = alpha,
                          gamma = gamma)

      if (!((opt$status >= 1) && (opt$status <= 4))) {
        warning("unsuccessful convergence, message ", opt$message,
                ", code ", opt$status, ", see documentation for nloptr")
      }

      optim.diagn <- list(convergence = opt$status, message = opt$message)

      list(optim.diagn = optim.diagn,
           theta = setNames(opt$theta, levels(regimes$regimes)),
           alpha = opt$alpha,
           sigma.sq = opt$sigma.sq,
           gamma = opt$gamma,
           loglik = opt$loglik)
    }
  )
)
