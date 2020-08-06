#' Significance Calculation of A Set of Generalized U Statistics.
#'
#'\code{zhangjtMat} calculates a matrix of p values for a set of generatlized U
#'statistics, with each element corresponding to a generalized U statistic.
#'
#' @param E the similarity matrix of one of the variates or multivariates.
#' @param G the similarity matrix of the other one of the variates or
#'  multivariates.
#' @param seq0 a set of powers for the similarity matrixes.
#'
#' @return a matrix of p values.
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
#' seq1 = c(0.5,1,2)
#' zhangjtMat(Kx,Ky,seq1)
zhangjtMat <- function(E,G,seq0){
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

  numKM = length(seq0)
  mat0 = matrix(NA,numKM,numKM)

  for(e11 in 1:numKM){
    for(g11 in 1:numKM){
      e1 = seq0[e11]
      g1 = seq0[g11]

      EVa1 = diag(EVa^e1)
      E1 = EV %*% EVa1 %*% t(EV)

      GVa1 = diag(GVa^g1)
      G1 = GV %*% GVa1 %*% t(GV)

      mat0[e11,g11] = zhangjtAppro(E1, G1)

    }
  }
  return(mat0)
}




#' Significance calculation for the combined generalized U statistc.
#'
#' @param mat0 a matrix of p values to be combined.
#'
#' @return the p value of the combined generalized U statistic.
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
#' seq1 = c(0.5,1,2)
#' mat1 = zhangjtMat(Kx,Ky,seq1)
#' CGU9(mat1)
CGU9 <- function(mat0){
  mat1 = tan((0.5-mat0)*pi)

  k = dim(mat0)[1]

  sta1 = sum(mat1)/(k^2)

  pval = 0.5 - atan(sta1)/pi

  return(pval)

}
