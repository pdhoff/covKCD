#' Kronecker-core covariance decomposition
#'
#' Computes the Kronecker-core decomposition of a covariance matrix. 
#'
#' The Kronecker-core decomposition is a representation of an 
#' arbitrary covariance matrix S in terms of a separable 
#' Kronecker covariance 
#' matrix K and a complementary non-separable core covariance matrix
#' C. The Kronecker covariance is the separable covariance matrix
#' that is closest to S
#' in terms of the divergence function
#' \deqn{ \log|K| + \mbox{trace}(K^{-1}S).} 
#' The core covariance matrix 
#' C is computed from 
#' S and K via 
#' \deqn{ C = K^{-1/2} S K^{-1/2}.}
#' 
#' @param S a covariance matrix of dimension (p1p2)*(p1p2).
#' @param p1 the row dimension.
#' @param p2 the column dimension. 
#' @param tol the convergence tolerance of the iterative algorithm. 
#' 
#' @return \code{covKCD} returns a list with the following elements: 
#' \describe{
#' \item{K}{the Kronecker covariance matrix;}
#' \item{C}{the core covariance matrix;}
#' \item{K1}{the row covariance matrix;}
#' \item{K2}{the column covariance matrix;} 
#' \item{div}{ the divergence between \code{S} and \code{K} across iterations of the algorithm.} 
#' }
#' @author Peter Hoff 
#' 
#' @examples 
#' p1<-4 ; p2<-3 ; n<-200
#' 
#' # create a matrix Y with separable covariance
#' A<-matrix(rnorm(p1*p1),p1,p1)  
#' B<-matrix(rnorm(p2*p2),p2,p2)/3
#' Y<-array(rnorm(n*p1*p2),dim=c(n,p1,p2))  
#' Y<-aperm( apply(Y,c(1,3),function(y){ A%*%y } ),c(2,1,3)) 
#' Y<-aperm( apply(Y,c(1,2),function(y){ B%*%y } ),c(2,3,1)) 
#' 
#' # covariance 
#' S<-mcov(Y) 
#' 
#' KCD<-covKCD(S,p1,p2) 
#' 
#' plot(A%*%t(A), KCD$K1)
#' plot(B%*%t(B), KCD$K2)
#' 
#' @export
covKCD<-function(S,p1,p2,tol=1e-8){  

  if( !( length(dim(S))==2 & all(dim(S)==p1*p2) ) ) {
    stop("covariance not correctly specified") 
  } 

  ## covariance array
  A<-covKCD::cm2ca(S,p1,p2) 
 
  ## starting values  
  K1<-matrix(0,p1,p1) ; K2<-matrix(0,p2,p2) 
  for(i in 1:p1){  K2<-K2+A[i,,i,]/p1 }
  for(j in 1:p2){  K1<-K1+A[,j,,j]/p2 }
  iK1<-solve(K1) ; iK2<-solve(K2) 

  ## starting objective function
  div1<-p1*log(det(K2)) + p2*log(det(K1)) + sum(S*kronecker(iK2,iK1))
  div0<-div1/2
  div<-NULL 

  ## iterative algorithm
  while( abs(div0-div1)/abs(div1) > tol ){

    ## update K1
    K1<-matrix(0,p1,p1)  
    for(i1 in 1:p1){for(i2 in i1:p1){ 
      K1[i1,i2]<-K1[i2,i1]<-sum(A[i1,,i2,]*iK2)/p2  }}
    iK1<-solve(K1) 

    ## update K2
    K2<-matrix(0,p2,p2)  
    for(j1 in 1:p2){for(j2 in j1:p2){ 
      K2[j1,j2]<-K2[j2,j1]<-sum(A[,j1,,j2]*iK1)/p1  }}
    iK2<-solve(K2)  

    ## compute divergence 
    div<-c(div,div1)
    div0<-div1 
    div1<-p1*log(det(K2)) + p2*log(det(K1)) + sum(S*kronecker(iK2,iK1))

  }

  ## compute core 
  Kmh<-kronecker(msqrt(iK2),msqrt(iK1)) 
  C<-Kmh%*%S%*%t(Kmh) 

list(K=kronecker(K2,K1),C=C,K1=K1,K2=K2,div=div) 
} 

