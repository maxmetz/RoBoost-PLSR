dummy = function (Y)
{
  if (sum(is.na(Y)) > 0)
    stop("NA in 'Y' are not allowed")
  Y <- as.factor(Y)
  lev <- levels(Y)
  nlev <- length(lev)
  if (nlev == 1)
    Y <- factor(Y, levels = c(lev, ".NA"))
  z <- model.matrix(~Y - 1)
  attr(z, "assign") <- NULL
  attr(z, "contrasts") <- NULL
  colnames(z) <- substring(colnames(z), first = 2, last = 1000000L)
  z
}



pls.nipalsw = function (X, Y, ncomp, weights = NULL)
{

  X <- .matrix(X)
  zdim <- dim(X)
  n <- zdim[1]
  zp <- zdim[2]
  Y <- .matrix(Y, row = FALSE, prefix.colnam = "y")
  q <- dim(Y)[2]
  if (is.null(weights))
    weights <- rep(1/n, n)
  else weights <- weights/sum(weights)

  xmeans <- NULL

  ymeans <- NULL


  nam <- paste("comp", 1:ncomp, sep = "")
  T <- matrix(nrow = n, ncol = ncomp, dimnames = list(row.names(X),
                                                      nam))
  R <- W <- P <- matrix(nrow = zp, ncol = ncomp, dimnames = list(colnames(X),
                                                                 nam))
  C <- matrix(nrow = q, ncol = ncomp, dimnames = list(colnames(Y),
                                                      nam))
  TT <- vector(length = ncomp)
  for (a in 1:ncomp) {
    XY <- crossprod(weights * X, Y)
    if (q == 1) {
      w <- XY
      w <- w/sqrt(sum(w * w))
    }
    else {
      w <- svd(XY, nu = 1, nv = 0)$u
    }
    t <- X %*% w
    tt <- sum(weights * t * t)
    c <- crossprod(weights * Y, t)/tt
    p <- crossprod(weights * X, t)/tt
    X <- X - tcrossprod(t, p)
    Y <- Y - tcrossprod(t, c)
    T[, a] <- t
    W[, a] <- w
    P[, a] <- p
    C[, a] <- c
    TT[a] <- tt
  }
  R <- W %*% solve(crossprod(P, W))
  list(T = T, P = P, W = W, C = C, R = R, TT = TT, xmeans = xmeans,
       ymeans = ymeans, weights = weights, T.ortho = TRUE,Y = Y,X = X)
}


.matrix = function (X, row = TRUE, prefix.colnam = "x")
{
  if (is.vector(X))
    if (row)
      X <- matrix(X, nrow = 1)
    else X <- matrix(X, ncol = 1)
    if (!is.matrix(X))
      X <- as.matrix(X)
    if (is.null(row.names(X)))
      row.names(X) <- 1:dim(X)[1]
    if (is.null(colnames(X)))
      colnames(X) <- paste(prefix.colnam, 1:dim(X)[2], sep = "")
    X
}

.center = function (X, center = matrixStats::colMeans2(X)) {
  t((t(X) - c(center)))}

.xmean = function (X, weights = NULL, row = FALSE)
{
  X <- .matrix(X, row = row)
  n <- dim(X)[1]
  if (is.null(weights))
    weights <- rep(1/n, n)
  else weights <- weights/sum(weights)
  colSums(weights * X)
}

.projscor = function (fm, X)
{
  T <- .center(.matrix(X), fm$xmeans) %*% fm$R
  rownam <- row.names(X)
  colnam <- paste("comp", 1:dim(T)[2], sep = "")
  dimnames(T) <- list(rownam, colnam)
  T
}

predict.pls = function(Xu,fm){
  Tu <- .projscor(fm, .matrix(Xu))
  m <- dim(Tu)[1]
  ncomp <- dim(Tu)[2]
  q = length(fm$ymeans)
  rownam.Xu <- row.names(Tu)
  Ymeans <- matrix(rep(fm$ymeans, m), nrow = m, byrow = TRUE)
  beta <- t(fm$C)
  for (a in 1:ncomp) {
    fit <- Ymeans + Tu[, 1:a, drop = FALSE] %*%
      beta[1:a, ,drop = FALSE]
  }
  return(list(fit = fit,Tu = Tu))
}


.findmax = function (x, seed = NULL)
{
  x <- which.max(x)
  n <- length(x)
  if (n > 1) {
    set.seed(seed = seed)
    x <- x[sample(seq_len(n), 1)]
    set.seed(seed = NULL)
  }
  x
}

F.weight = function(x,cw = 4){
  if(cw == Inf){x = rep(1,length(as.vector(x)))}else{
    e <-as.vector(x)

    s = median(abs(x))

    x = e/(cw*s)
    index = c(1:length(x))[abs(x)>=1]
    if(length(index)==0){ x[] = 1
    }else{
      x[-c(index)] = (1-x[-c(index)]^2)^2
      x[c(index)] = 0}}

  return(x)
}



rmse.t = function(x, W ,z){
  o = order(W)
  p = z*length(o)
  o = o[1:p]
  x = x[-o]
  msep <- sum((x^2)/length(x))
  rmsep <- sqrt(msep)
  return(rmsep)}

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


.estim_poid = function (x, cw, med)
{
  if (cw == Inf) {
    x = rep(1, length(as.vector(x)))
  }
  else {
    e <- as.vector(x)
    s = med
    x = e/(cw * s)
    index = c(1:length(x))[abs(x) >= 1]
    if (length(index) == 0) {
      x[] = 1
    }
    else {
      x[-c(index)] = (1 - x[-c(index)]^2)^2
      x[c(index)] = 0
    }
  }
  return(x)
}


.predict_rob_outlier = function(mod, Xu, Yu)
{
  if (!is(mod, "roboost_PLSR")) {
    stop("This is not an RoBoost-PLSR object")
  }
  ncomp = length(mod$fm)
  m = nrow(Xu)
  q = nrow(mod$C)
  Ymeans <- matrix(rep(mod$fm[[1]]$ymeans, m), nrow = m, byrow = TRUE)
  pred <- array(dim = c(m, ncomp + 1, q))
  pred[, 1, ] <- Ymeans
  beta.p <- t(mod$C)
  Tu = .center(.matrix(Xu), mod$fm[[1]]$xmeans) %*% mod$R
  for (a in seq_len(ncomp)) {
    pred[, a + 1, ] <- Ymeans + Tu[, seq_len(a), drop = FALSE] %*%
      beta.p[seq_len(a), , drop = FALSE]
  }

  nam <- paste("comp", 1:ncomp, sep = "")
  a = matrix(nrow = nrow(Xu), ncol = ncomp, dimnames = list(row.names(X), nam))
  b = matrix(nrow = nrow(Xu), ncol = ncomp, dimnames = list(row.names(X), nam))
  g = matrix(nrow = nrow(Xu), ncol = ncomp, dimnames = list(row.names(X), nam))

  fit. = c()
  for ( i in 1:ncomp) {
    r = Yu - pred[,i+1,]
    r1 = .center(Xu, mod$fm[[i]]$xmeans) - Tu[,1:i, drop=FALSE] %*% t(mod$P[,1:i,drop = FALSE])
    r1 = sqrt(r1*r1)
    r1 = rowSums(r1)
    l = Tu[,i]

    a[,i] = .estim_poid(r1, mod$hyp_para$alpha, mod$crit$r1[i])
    b[,i] = .estim_poid(r, mod$hyp_para$beta, mod$crit$r[i])
    g[,i] = .estim_poid(l, mod$hyp_para$gamma, mod$crit$l[i])

    fit = as.data.frame(r)
    fit$ncomp = rep(i,nrow(fit))
    fit$Yt = Yu
    fit$b = .estim_poid(r1, mod$hyp_para$alpha, mod$crit$r1[i])*.estim_poid(r, mod$hyp_para$beta, mod$crit$r[i])*.estim_poid(l, mod$hyp_para$gamma, mod$crit$l[i])
    fit. = rbind(fit.,fit)
  }
  fit.$ncomp = as.factor(fit.$ncomp)
  W_est = list(Wx = a, Wy = b, PL = g)
  fm = list(W_est = W_est, fit = fit., pred = pred, Tu = Tu)
  return(fm)
}

