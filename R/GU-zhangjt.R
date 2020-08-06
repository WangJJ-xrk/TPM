#' Similarity Matrix Calculation.
#'
#' \code{calcKerMat_Cen} calculates a centered similarity matrix for a data
#' matrix.
#'
#' @param X the data matrix, each row represents an observation.
#'
#' @return the centered similarity matrix.
#'
#' @importFrom stats dist
#'
#' @export
#'
#' @examples
#' library(MASS)
#' n = 50
#' p = 100
#' k = 15
#' sigmax = diag(rep(0.5,k)) + matrix(0.5,k,k)
#' x = mvrnorm(n, rep(0,k), sigmax)
#' Kx = calcKerMat_Cen(x)
calcKerMat_Cen <- function(X){

  dis.X = as.matrix(dist(X))

  K.X = -0.5*dis.X^2

  n = dim(K.X)[1]

  H = diag(rep(1,n)) - matrix(1/n,n,n)

  cK.X = H %*% K.X %*% H

  return(cK.X)

}


#' Significance Calculation for the Chi-squared-type Mixtures.
#'
#'\code{zhangjtAppro} calculates the p value of the chi-squared-type mixtures
#'via the method of Jinting Zhang (2013).
#'
#' @param E the similarity matrix of one of the variates or multivariates.
#' @param G the similarity matrix of the other one of the variates or
#'  multivariates.
#'
#' @return the p value of the generalized U statistic, i.e.,
#'  $\frac 1n \text {tr}(EG)$, where $n$ is the sample size.
#' @export
#'
#' @examples
#' library(MASS)
#' n = 50
#' p = 100
#' k = 15
#' sigmax = diag(rep(0.5,k)) + matrix(0.5,k,k)
#' sigmay = diag(rep(1,p))
#' for(i in 1:p){
#'   for(j in 1:p){
#'     sigmay[i,j] = 0.5^abs(i-j)
#'   }
#' }
#' r1 = 0.05
#' beta0 = r1*matrix(rbinom(k*p,1,0.9), k, p)
#' x = mvrnorm(n, rep(0,k), sigmax)
#' y = x%*%beta0 + mvrnorm(n, rep(0,p), sigmay)
#' Kx = calcKerMat_Cen(x)
#' Ky = calcKerMat_Cen(y)
#' zhangjtAppro(Kx,Ky)
zhangjtAppro <- function(E,G){
  n = dim(E)[1]

  eigE = eigen(E)
  Eval = eigE$values
  Evec = eigE$vectors
  kqe = sum(Eval > 0.0001)
  EV = Evec[,1:kqe]
  EVa = Eval[1:kqe]

  eigG = eigen(G)
  Gval = eigG$values
  Gvec = eigG$vectors
  kq = sum(Gval > 0.0001)
  GV = Gvec[,1:kq]
  GVa = Gval[1:kq]

  Cr0 = outer(EVa, GVa, "*")
  Cr = Cr0/(n^2)

  s1 = sum(Cr)
  s2 = sum(Cr^2)
  s3 = sum(Cr^3)

  alpha0 = s3/s2
  beta0 = s1 - s2^2/s3
  d0 = s2^3/s3^2

  T1 = sum(E*G)/n - beta0
  p1 = 1 - pgamma(T1, shape = d0/2, scale = 2*alpha0)
  return(p1)
}

