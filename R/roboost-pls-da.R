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
roboost_plsr <- function (X,Y,ncomp,niter = 50,gamma =gamma,beta=beta,alpha=alpha,th = 1-10^-12){

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


  a = 1
  fm = list()
  Y1 = Y
  X1 = X


  while(a<(ncomp+1)){
    cor = 0
    th = th
    f = 1
    if(a==1){d =  rep(1, nrow(X))}


    while(cor<th){

      if(a==1){
        xmeans <- .xmean(X1, weights = d)
        X <- .center(X1, xmeans)
        ymeans <- .xmean(Y1, weights = d)
        Y <- .center(Y1, ymeans)
      }else{
        xmeans = rep(0,ncol(X))
        ymeans = rep(0,ncol(Y))
      }



      fm[[a]]  = pls.nipalsw(X,Y,ncomp = 1,weights = d)
      fm[[a]]$xmeans = xmeans
      fm[[a]]$ymeans = ymeans


      r = fm[[a]]$Y
      #plot(r[,1])


      for (i in 1:ncol(r)){
        r[,i]  = r[,i]
        r[,i] = F.weight(x = r[,i],cw = beta)

      }

      #plot(r[,1])
      #plot(r[,2])
      beta.w = r = apply(r,1,FUN = prod)
      r = r/sum(r)
      ##plot(r)

      r1 = (X - tcrossprod(fm[[a]]$T, fm[[a]]$P))
      r1 = sqrt(rowSums(r1 * r1))
      r1  = r1

      alpha.w = r1 <- F.weight(r1,alpha)
      r1 = as.vector(r1)
      r1 = r1/sum(r1)
      #plot(r1)

      l = fm[[a]]$T
      l = l
      gamma.w = l = F.weight(l,gamma)
      l = l/sum(l)
      #plot(l)

      d <- fm[[a]]$weights/sum(fm[[a]]$weights)
      d <- (r*r1*l)/sum(r*r1*l)
      #plot(d)

      u = as.matrix(Y)%*%as.matrix(fm[[a]]$C)
      q = crossprod(d * fm[[a]]$T, u)/ sum(fm[[a]]$T * u * fm[[a]]$T)
      #print(q)

      if(f>1){
        cor = min(c(q,q1))/max(c(q,q1))}

      if(f>niter) {cor = Inf}
      f = f+1
      q1 = q

    }

    w = fm[[a]]$weights
    fm[[a]]$list.w = list(alpha.w,beta.w,gamma.w)
    X = fm[[a]]$X
    Y = fm[[a]]$Y

    T[, a] <- fm[[a]]$T
    W[, a] <- fm[[a]]$W
    P[, a] <- fm[[a]]$P
    C[, a] <- fm[[a]]$C
    a = a+1
  }
  R =  W %*% solve(crossprod(P, W))
  class(fm) = "roboost_PLSR"
  fm = list(T = T, P = P, W = W, C = C, R = R, fm = fm)
  class(fm) = "roboost_PLSR"
  return(fm)

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
predict_roboost_plsr = function(mod,Xu){
  if(!is(mod,"roboost_PLSR"))
  {stop("This is not an RoBoost-PLSR object")}


  ncomp = length(mod)
  fit.tot = list()
  for( i in (1:ncomp)){
    pred = predict.pls(mod[[i]],Xu = Xu)
    X.hat =   pred$Tu%*%t(mod[[i]]$P)
    Xu = .center(Xu,mod[[i]]$xmeans)
    Xu = Xu - X.hat
    if(i!=1) {fit.tot[[i]] = fit.tot[[i-1]] + pred$fit
    }else{fit.tot[[i]] = pred$fit}

  }
  return(fit.tot)

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
predict_roboost_plsr_rot = function(mod,Xu){
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
roboost_bcoef = function(mod,ncomp){

  beta <- t(mod$C)[seq_len(ncomp), , drop = FALSE]
  b <- mod$R[, seq_len(ncomp), drop = FALSE] %*% beta
  int <- mod$fm[[1]]$ymeans - t(mod$fm[[1]]$xmeans) %*% b
  b <- rbind(int, b)
  row.names(b)[1] <- "intercept"
  b
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


roboost_plda <- function (X,Y,ncomp,niter = 50,gamma =gamma,beta=beta,alpha=alpha,th = 1-10^-12,init.weights=NULL){

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


  a = 1
  fm = list()
  Y1 = Y
  X1 = X
  Y.logical = apply(Y,2,as.logical)

  while(a<(ncomp+1)){
    cor = 0
    th = th
    f = 1
    if(a==1&is.null(init.weights)){d =  rep(1, nrow(X))}
    if(a==1&!is.null(init.weights)){d = init.weights;print("cest ok")}

    while(cor<th){

      if(a==1){
        xmeans <- .xmean(X1, weights = d)
        X <- .center(X1, xmeans)
        ymeans <- .xmean(Y1, weights = d)
        Y <- .center(Y1, ymeans)
      }else{
        xmeans = rep(0,ncol(X))
        ymeans = rep(0,ncol(Y))
      }



      fm[[a]]  = pls.nipalsw(X,Y,ncomp = 1,weights = d)
      fm[[a]]$xmeans = xmeans
      fm[[a]]$ymeans = ymeans


      fit = Y1 - fm[[a]]$Y
      fit = exp(fit)/rowSums(exp(fit))

      r = Y1 - fit

      d <- fm[[a]]$weights/sum(fm[[a]]$weights)

      for (i. in 1:ncol(r)){
        for (i in 1:ncol(Y.logical)){

          x = r[c(Y.logical[,i]),i.]- median(r[c(Y.logical[,i]),i.])
          r[c(Y.logical[,i]),i.] = F.weight(x,beta)

        }
      }


      beta.w = r = apply(r,1,FUN = prod)
      r = r/sum(r)

      r1 = fm[[a]]$X
      r1 = as.vector(sqrt(rowSums(r1 * r1)))

      for (i in 1:ncol(Y.logical)){
        x = r1[c(Y.logical[,i])]-median(r1[c(Y.logical[,i])])
        r1[c(Y.logical[,i])] = F.weight(x,alpha)
      }
      alpha.w = r1 = r1/sum(r1)


      l = vector(length = nrow(Y))
      for (i in 1:ncol(Y.logical)){

        x = fm[[a]]$T[c(Y.logical[,i]),]-median(fm[[a]]$T[c(Y.logical[,i]),])
        l[c(Y.logical[,i])] = F.weight(x,gamma)

      }

      gamma.w = l = l/sum(l)

      d <- fm[[a]]$weights/sum(fm[[a]]$weights)
      d <- (r1*r*l)/sum(r1*r*l)

      u = as.matrix(Y)%*%as.matrix(fm[[a]]$C)
      q = crossprod(d * fm[[a]]$T, u)/ sum(fm[[a]]$T * u * fm[[a]]$T)


      if(f>1){
        cor = min(c(q,q1))/max(c(q,q1))}

      if(f>niter) {cor = Inf}
      f = f+1
      q1 = q

    }

    w = fm[[a]]$weights
    fm[[a]]$list.w = list(alpha.w,beta.w,gamma.w)
    X = fm[[a]]$X
    Y = fm[[a]]$Y

    T[, a] <- fm[[a]]$T
    W[, a] <- fm[[a]]$W
    P[, a] <- fm[[a]]$P
    C[, a] <- fm[[a]]$C
    a = a+1
  }
  R =  W %*% solve(crossprod(P, W))
  class(fm) = "roboost_PLSR"
  fm = list(T = T, P = P, W = W, C = C, R = R, fm = fm)
  class(fm) = "roboost_PLSR"
  return(fm)

}




