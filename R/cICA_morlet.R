#'
#'
#' This function performs colorICA Morlet algorithm on CIFTI files.
#' @import polspline
#' @import parallel
#' @import dplR
#'
#'
#' @param Xin CIFTI file data
#' @param M Xin first dimension
#' @param Win diagonal of matrix M
#' @param tol tolerance
#' @param maxit maximum number of iteration
#' @param nmaxit maximum number of iteration
#' @param dj Divide octave in sub-octaves. If dj = 0.25 this will do 4 sub-octaves per octave.
#' @param ncores Number of cores to use during fitting
#' @return  return(as.list(result))

#' @export

cICA_morlet<-function (Xin,
                   M = dim(Xin)[1],
                   Win = diag(M),
                   tol = 1e-04,
                   maxit = 20,
                   nmaxit=5,
                   dj=0.25,
                   ncores=2)

  #sopoc : Sub-octaves per octave calculated.
{
  p = dim(Xin)[1]
  if (M > p) {
    stop("Number of sources must be less or equal than number \n  of variables")
  }

  ### pre-whitening
  T = ncol(Xin)
  Xc = t(scale(t(Xin), center = TRUE, scale = FALSE))
  svdcovmat = svd(Xc/sqrt(T))
  K = t(svdcovmat$u %*% diag(1/svdcovmat$d))
  K = K[1:M, ]
  Xc = K %*% Xc

 ### Initialization
  Wold = Wnew = Wtmp = Win
  X.morlet = lapply(1:nrow(Xc),function(idx){morlet(y1 = Xc[idx,], x1 = 1:T, dj = dj, siglvl = 0.95)})
  perid.x = X.morlet[[1]]$period
  S.mat = Wold %*% Xc
  S.mat.ppt = mvfft(t(S.mat))
  freqlength=floor(T/2)
  nperid = length(perid.x)
  fft.indx=2:(freqlength)

  lspec.ests = apply(S.mat.ppt, 2,function(obj){
        mm = length(obj)
        tmpfit=lspec(period=Mod(fft(obj)[2:(mm/2)])^2/(2*pi*mm))
        pred = dlspec( 1/perid.x *2*pi,tmpfit)
        lspec.est = pred$d + pred$m ## need to straightedn up.
        return(lspec.est)})
  MtX = do.call(cbind,lapply(1:length(perid.x),
               function(k){
                 tmp1 = do.call(cbind, lapply(X.morlet, function(x)Re(x$wave[,k])))
                 tmp2 = do.call(cbind, lapply(X.morlet, function(x)Im(x$wave[,k])))
                 re=as.vector( t(tmp1) %*% Conj(tmp1)/T + t(tmp2) %*% Conj(tmp2)/T)
                 return(re)}))

  Ixf = MtX %*% (1/lspec.ests)
  ### Optimization
  lim=1
  iter=0
  wlik = -Inf
  NInv = 0
  Ndec = 0

  while (lim > tol & iter < maxit & NInv < nmaxit) {
    iter = iter + 1
    taucount = 1
    err = 1
    orthoerror = 1
    tau = 0.5
    eigenval = rep(0, M)
    Wold = Wnew
    ## update W
    while (taucount < 60 & err > 1e-05 & orthoerror > 1e-05) {
      for (j in 1:M) {
        Gam = 0
        if (j > 1) {
          for (k in 1:(j - 1)) {
            nu = matrix(Wnew[k, ], M, 1) %*% matrix(Wnew[k,], 1, M)
            Gam = Gam + nu
          }
        }

        tmpV = matrix(Ixf[,j],M,M) + tau * Gam
        eigenv = eigen(tmpV)
        eigenval[j] = eigenv$values[M]
        Wnew[j, ] = eigenv$vectors[,M]
      }
      orthoerror = (sum(sum((Wnew %*% t(Wnew) - diag(rep(1, M)))^2)))
      err = amari_distance(rerow(Wnew), rerow(Wold))
      taucount = taucount + 1
      tau = 2 * tau
    }

    ## Update spectral density (when iter=1, it is initialization)
    S.mat = Wnew %*% Xc
    S.mat.ppt = mvfft(t(S.mat))
    freqlength=floor(T/2)
    nperid = length(perid.x)
    fft.indx=2:(freqlength)

    lspec.ests = apply(S.mat.ppt, 2,function(obj){
      mm = length(obj)
      tmpfit=lspec(period=Mod(fft(obj)[2:(mm/2)])^2/(2*pi*mm))
      pred = dlspec( 1/perid.x *2*pi,tmpfit)
      lspec.est = pred$d + pred$m ## need to straightedn up.
      return(lspec.est)})

    MtX = do.call(cbind,lapply(1:length(perid.x),
                               function(k){
                                 tmp1 = do.call(cbind, lapply(X.morlet, function(x)Re(x$wave[,k])))
                                 tmp2 = do.call(cbind, lapply(X.morlet, function(x)Im(x$wave[,k])))
                                 re=as.vector( t(tmp1) %*% Conj(tmp1)/T + t(tmp2) %*% Conj(tmp2)/T)
                                 return(re)}))

    Ixf = MtX %*% (1/lspec.ests)

    ## Whittle likelihood
    wlik2 = -1 * sum(eigenval) - 1 * do.call(sum,lapply(lspec.ests, function(x)sum(log(x)))) + T*log(abs(det(Wnew)))

    lim = err
    print(paste("cICA-tf - Iteration ", iter, ": error is equal to ",lim, sep = ""))


    if (wlik < wlik2 ) {
      #       print('Whittle likelihood increased.')
      wlik = wlik2
      Wtmp = Wnew
    }else{print(paste("cICA-tf - Iteration ", iter,
                      ": current Whittle likelihood(", wlik2, ") is smaller than previous one (",
                      wlik, ")."))
      Wtmp = Wold ## Wtmp is the temporal saving of the local maxima
      Ndec = Ndec+1 ## we will allow 3 consecutive decrease in the whittle likelihood.
    }

    if (( (iter == maxit |  Ndec==4) & NInv < nmaxit)) {
      print("Color ICA: iteration reaches to maximum. Start new iteration.")
      Wnew = qr.Q(qr(matrix(rnorm(M * M), M, M)))
      iter = 0
      Ndec = 0
      NInv = NInv + 1
      lim=1
    }

    if (NInv == nmaxit) {
      print("Color ICA: no convergence")
    }

    ## update old W to new
  }

  wt = Wtmp %*% K
  result = new.env()
  result$W = Wtmp
  result$K = K
  result$A = t(wt) %*% solve(wt %*% t(wt))
  result$S = wt %*% Xin
  result$X = Xin
  result$iter = iter
  result$NInv = NInv
  as.list(result)
}



