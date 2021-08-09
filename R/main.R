#' RoBoost-PLSR : Robust method for partial least squares regression
#'
#' @param X
#' Explanatory variables
#' @param Y
#' Explained Variables
#' @param ncomp
#' Number of latent variables
#' @param niter
#' Number of maximal iterations
#' @param gamma
#' parameters for leverage point
#' @param beta
#' parameters for Y-residuals
#' @param alpha
#' parameters for X-residuals
#'
#'
#' @return
#' @export
#'
#' @examples
#'
#' n <- 10
#' p <- 6
#' set.seed(1)
#' X <- matrix(rnorm(n * p, mean = 10), ncol = p)
#' y1 <- 100 * rnorm(n)
#' y2 <- 100 * rnorm(n)
#' Y <- cbind(y1, y2)
#' set.seed(NULL)

#' Xr <- X[1:8, ] ; Yr <- Y[1:8, ]
#' Xu <- X[9:10, ] ; Yu <- Y[9:10, ]


#' library(roboost)

#' ncomp = 3

#' alpha = Inf
#' beta  = Inf
#' gamma = Inf

#' mod = roboost_plsr(X = Xr,Y = Yr ,ncomp,gamma =gamma,beta=beta,alpha=alpha)
#' pred = predict_roboost_plsr(mod$fm,Xu)
RoBoost_plsr <- function (X, Y, ncomp, niter = 50, gamma = gamma, beta = beta,
                          alpha = alpha, th = 1 - 10^-12)
{
  zdim <- dim(X)
  n <- zdim[1]
  zp <- zdim[2]
  Y <- .matrix(Y, row = FALSE, prefix.colnam = "y")
  q <- dim(Y)[2]
  xmeans <- NULL
  ymeans <- NULL
  nam <- paste("comp", 1:ncomp, sep = "")
  T <- matrix(nrow = n, ncol = ncomp, dimnames = list(row.names(X),
                                                      nam))
  R <- W <- P <- matrix(nrow = zp, ncol = ncomp, dimnames = list(colnames(X),
                                                                 nam))
  C <- matrix(nrow = q, ncol = ncomp, dimnames = list(colnames(Y),
                                                      nam))
  crit = list()
  Cr = matrix(nrow = 1, ncol = ncomp, dimnames = list(row.names(X),
                                                      nam))
  Cr1 = matrix(nrow = 1, ncol = ncomp, dimnames = list(row.names(X),
                                                       nam))
  Cl = matrix(nrow = 1, ncol = ncomp, dimnames = list(row.names(X),
                                                      nam))
  a = 1
  fm = list()
  Y1 = Y
  X1 = X
  while (a < (ncomp + 1)) {
    cor = 0
    th = th
    f = 1
    if (a == 1) {
      d = rep(1, nrow(X))
    }
    while (cor < th) {
      if (a == 1) {
        xmeans <- .xmean(X1, weights = d)
        X <- .center(X1, xmeans)
        ymeans <- .xmean(Y1, weights = d)
        Y <- .center(Y1, ymeans)
      }
      else {
        xmeans = rep(0, ncol(X))
        ymeans = rep(0, ncol(Y))
      }
      #if(f>1){
      #cr = median(abs(ry))
      #cr1 = median(abs(rx))
      #cl = median(abs(rl))
      #beta.w = F.weight(ry, beta)
      #alpha.w = F.weight(rx, alpha)
      #gamma.w = F.weight(rl, gamma)
      #}

      fm[[a]] = pls.nipalsw(X, Y, ncomp = 1, weights = d)
      fm[[a]]$xmeans = xmeans
      fm[[a]]$ymeans = ymeans
      ry = r = fm[[a]]$Y
      for (i in 1:ncol(r)) {
        r[, i] = r[, i]
        r[, i] = F.weight(x = r[, i], cw = beta)
      }
      r = apply(r, 1, FUN = prod)
      r = r/sum(r)
      r1 = (X - tcrossprod(fm[[a]]$T, fm[[a]]$P))
      r1 = sqrt(rowSums(r1 * r1))
      rx = r1
      r1 <- F.weight(r1, alpha)
      r1 = as.vector(r1)
      r1 = r1/sum(r1)
      l = fm[[a]]$T
      rl = l
      l = F.weight(l, gamma)
      l = l/sum(l)
      d <- fm[[a]]$weights/sum(fm[[a]]$weights)
      d <- (r * r1 * l)/sum(r * r1 * l)
      u = as.matrix(Y) %*% as.matrix(fm[[a]]$C)
      q = crossprod(d * fm[[a]]$T, u)/sum(fm[[a]]$T *
                                            u * fm[[a]]$T)
      if (f > 1) {
        cor = min(c(q, q1))/max(c(q, q1))
      }
      if (f > niter) {
        cor = Inf
      }
      f = f + 1
      q1 = q
    }
    beta.w = F.weight(ry, beta)
    alpha.w = F.weight(rx, alpha)
    gamma.w = F.weight(rl, gamma)
    w = fm[[a]]$weights
    fm[[a]]$list.w = list(alpha.w, beta.w, gamma.w)
    X = fm[[a]]$X
    Y = fm[[a]]$Y
    T[, a] <- fm[[a]]$T
    W[, a] <- fm[[a]]$W
    P[, a] <- fm[[a]]$P
    C[, a] <- fm[[a]]$C
    Cr[,a] <- median(abs(ry))
    Cr1[,a] <- median(abs(rx))
    Cl[,a] <- median(abs(rl))
    a = a + 1
  }
  hyppara = list(alpha= alpha, beta = beta, gamma = gamma)
  crit$r = Cr
  crit$r1 = Cr1
  crit$l = Cl
  R = W %*% solve(crossprod(P, W))
  class(fm) = "roboost_PLSR"
  fm = list(T = T, P = P, W = W, C = C, R = R, fm = fm, crit = crit, hyp_para = hyppara)
  class(fm) = "roboost_PLSR"
  return(fm)
}

#' RoBoost-PLSR : Robust method for partial least squares regression
#'
#' @param X
#' Explanatory variables
#' @param Y
#' Explained Variables
#' @param ncomp
#' Number of latent variables
#' @param niter
#' Number of maximal iterations
#' @param gamma
#' parameters for leverage point
#' @param beta
#' parameters for Y-residuals
#' @param alpha
#' parameters for X-residuals
#' @param th
#' convergence threshold
#' @param segm
#' list of index
#'
#' @return
#' @export
#' @examples
#'
#'
#' n <- 100
#' p <- 10
#' set.seed(1)
#' X <- matrix(rnorm(n * p, mean = 10), ncol = p)
#' Y <- 100 * rnorm(n)
#'
#' library(dplyr)
#' n =nrow(X)
#' segm = segmkf(n = n,K = 5, nrep = 5)
#' ncomp = 5
#' th = 1-10^-12
#' niter = 20
#' alpha = Inf
#' beta = 4
#' gamma = Inf
#' test = CVROB(X = X, Y = Y ,segm = segm,alpha = alpha,beta = beta,gamma = gamma,ncomp = ncomp,niter = niter,th = th)
#' rmsecv = test %>% group_by(ncomp) %>% summarise(rmse.t(r,b,z = 0.15))
#' colnames(rmsecv)[2] = "r"
#' rmsecv[rmsecv$r==min(rmsecv$r),]
#' r2cv = test %>% group_by(ncomp) %>% summarise(r2.t(r,Yt,b,z = 0.15))
#' colnames(r2cv)[2] = "r"
#' r2cv[r2cv$r==max(r2cv$r),]




CVROB = function(X,Y,segm, alpha, beta, gamma, ncomp, niter, th) {


  para = expand.grid(alpha,beta,gamma)
  for (j in 1 : nrow(para)) {

    p = para[j,]
    fit. = c()
    param = c()
    rep = c()


    for (i in 1:length(segm)) {
      for (t in 1:length(segm[[i]])) {

        mod = RoBoost_plsr(X[-segm[[i]][[t]],],Y[-segm[[i]][[t]]], ncomp = ncomp,niter = niter,alpha = p[[1]], beta = p[[2]], gamma = p[[3]], th = th)
        prediction = .predict_rob_outlier(mod,Xu = as.matrix(X[segm[[i]][[t]],drop=FALSE,]), Yu = Y[segm[[i]][[t]]])
        fit. = rbind(fit.,prediction$fit)
        rep = c(rep,rep(i,nrow(prediction$fit)))
        param = c(param,rep(j,nrow(prediction$fit)))

      }
    }
  }
  return(data.frame(fit., rep = rep, param = param))
}

#' Segments for cross-validation
#'
#'
#' @param
#' n : number of samples
#' @param
#' K : number of fold
#' @param
#' y : A vector of length n defining the blocks.
#' @param
#' rep : The number of replications of the repeated CV. Default to nrep = 1.
#' @param
#' For segmkf. The type K-fold CV. Possible values are "random" (default), "consecutive" and "interleaved".
#' @return
#' @export
#'
#' @examples
#'
#'

segmkf = function (n, y = NULL, K = 5, type = c("random", "consecutive",
                                                "interleaved"), nrep = 1, seed = NULL)
{
  type <- match.arg(type)
  segm <- vector("list", length = nrep)
  names(segm) <- paste("rep", seq_len(nrep), sep = "")
  n <- round(n)
  zn <- n
  if (!is.null(y)) {
    if (length(y) != n)
      stop("y must be of size n")
    yagg <- unique(y)
    zn <- length(yagg)
  }
  lseg <- ceiling(zn/K)
  nna <- K * lseg - zn
  set.seed(seed = seed)
  for (i in seq_len(nrep)) {
    z <- switch(type, random = matrix(c(sample(seq_len(zn)),
                                        rep(NA, nna)), ncol = K, byrow = TRUE), consecutive = {
                                          x <- c(matrix(c(rep(1, zn), rep(NA, nna)), ncol = K,
                                                        byrow = TRUE))
                                          x[!is.na(x)] <- cumsum(na.omit(x))
                                          x <- matrix(x, ncol = K, byrow = FALSE)
                                          x
                                        }, interleaved = matrix(c(seq_len(zn), rep(NA, nna)),
                                                                ncol = K, byrow = TRUE))
    z <- lapply(data.frame(z), FUN = function(x) c(na.omit(x)))
    names(z) <- paste("segm", seq_len(K), sep = "")
    segm[[i]] <- z
  }
  if (!is.null(y)) {
    vecn <- seq_len(n)
    zsegm <- segm
    for (i in seq_len(nrep)) {
      for (j in seq_len(K)) {
        u <- segm[[i]][[j]]
        v <- which(y %in% yagg[u])
        zsegm[[i]][[j]] <- v
      }
    }
    segm <- zsegm
  }
  set.seed(seed = NULL)
  segm
}

#' Prediction with RoBoost-PLSR
#'
#'
#' @param
#' mod : RoBoost-PLSR model
#' @param
#' Xu : New dataset
#'
#' @return
#' @export
#'
#' @examples
#'
#'
predict_RoBoost_plsr = function(mod,Xu){
  if(!is(mod,"roboost_PLSR"))
  {stop("This is not an RoBoost-PLSR object")}
  ncomp = length(mod$fm)
  fit.tot = list()
  m = nrow(Xu)
  q = nrow(mod$C)

  Ymeans <- matrix(rep(mod$fm[[1]]$ymeans, m), nrow = m, byrow = TRUE)
  fit <- array(dim = c(m, ncomp + 1, q))
  fit[, 1, ] <- Ymeans
  beta <- t(mod$C)
  Tu = .center(.matrix(Xu), mod$fm[[1]]$xmeans) %*% mod$R
  for (a in seq_len(ncomp)) {
    fit[, a + 1, ] <- Ymeans + Tu[, seq_len(a), drop = FALSE] %*%
      beta[seq_len(a), , drop = FALSE]
  }
  return(fit)

}


#' rmse trimmed
#'
#' @param x
#' residuals
#'
#' @param W
#' weigths
#' @param z
#' trimmed threshold
#' @return
#' @export
#'
#' @examples
#'
#'


rmse.t = function(x, W ,z){
  o = order(W)
  p = z*length(o)
  o = o[1:p]
  x = x[-o]
  msep <- sum((x^2)/length(x))
  rmsep <- sqrt(msep)
  return(rmsep)}


#' rmse trimmed
#'
#' @param x
#' residuals
#' @param y
#' Y responses
#' @param W
#' weigths
#' @param z
#' trimmed threshold
#' @return
#' @export
#'
#' @examples
#'
#'

r2.t = function(x,y,W,z){
  o = order(W)
  p = z*length(o)
  o = o[1:p]
  x = x[-o]
  y = y[-o]
  msep <- sum((x^2)/length(x))
  .foo <- function(x) {
    sum(((x - mean(x))^2))
  }
  y.msep = .foo(y)/length(x)
  r2 = (msep/y.msep)
  r2 = 1 - r2
  return(r2)}



#' bcoef for RoBoost-PLSR
#'
#'
#' @param
#' mod : RoBoost-PLSR model
#' @param
#' ncomp : the number of lv to estimate the model
#' @return
#' @export
#'
#' @examples
#'
#'
RoBoost_bcoef = function(mod,ncomp){

  beta <- t(mod$C)[seq_len(ncomp), , drop = FALSE]
  b <- mod$R[, seq_len(ncomp), drop = FALSE] %*% beta
  int <- mod$fm[[1]]$ymeans - t(mod$fm[[1]]$xmeans) %*% b
  b <- rbind(int, b)
  row.names(b)[1] <- "intercept"
  b
}
