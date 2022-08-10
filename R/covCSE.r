#' Log multivariate gamma function
#' 
#' Compute the logarithm of the multivariate gamma function 
#' \eqn{\log \Gamma_p(a)}. 
#' 
#' @param a a numeric scalar.
#' @param p a positive integer.
#' 
#' @return a scalar 
#' 
#' @author Peter Hoff
#' @export
lmvgamma<-function(a,p){  
   if( (as.integer(p)-p!=0) | p<1 ){ stop("p must be a positive integer") }
   .25*p*(p-1)*log(pi) + sum(lgamma( a+(1-(1:p))/2 ))
} 


#' Empirical Bayes core shrinkage covariance estimator
#'
#' Estimate a covariance matrix by adaptively shrinking the core.  
#' 
#' @param data either a numeric n*p1*p2 array consisting of n data matrices 
#' each of dimension p1*p2, or a p1*p2 covariance matrix of data of this type.
#' If the latter, the values of n, p1 and p2 must be specified.
#' @param n the sample size.
#' @param p1 the row dimension of the data matrices.
#' @param p2 the column dimension of the data matrices. 
#' @param tol the convergence tolerance of the iterative algorithm.  
#' 
#' @return a covariance matrix of the same dimension as \code{S}. 
#' The attribute \code{w} of \code{S} gives the shrinkage weight on the 
#' Kronecker covariance of \code{S}. 
#' 
#' @author Peter Hoff 
#'
#' @examples 
#' p1<-4 ; p2<-3 ; n<-20
#' 
#' # create a matrix Y with separable covariance
#' Sig1<-rWishart(1,p1,diag(p1))[,,1] 
#' Sig2<-rWishart(1,p2,diag(p2))[,,1] 
#' 
#' Y<-array(rnorm(n*p1*p2),dim=c(n,p1,p2))  
#' Y<-aperm( apply(Y,c(1,3),function(y){ msqrt(Sig1)%*%y } ),c(2,1,3)) 
#' Y<-aperm( apply(Y,c(1,2),function(y){ msqrt(Sig2)%*%y } ),c(2,3,1)) 
#' 
#' # covariance 
#' S<-mcov(Y) 
#' covCSE(S,n,p1,p2)  
#'
#' # now an unstructured covariance
#' S<-rWishart(1,p1*p2,diag(p1*p2))[,,1] 
#' covCSE(S,n,p1,p2) 
#' 
#' @import stats
#' @export
covCSE<-function(data,n=NULL,p1=NULL,p2=NULL,tol=1e-08){ 

  if(length(dim(data))==2 & length(c(n,p1,p2))==3 && all(dim(data)==p1*p2)){
    S<-data 
  } else 
  if(length(dim(data))==3){ 
    S<-covKCD::mcov(data) 
    n<-dim(data)[1] 
    p1<-dim(data)[2] 
    p2<-dim(data)[3]
  } else stop("data not correctly specified") 

  # Kronecker-core decomposition
  KCD<-covKCD::covKCD(S,p1,p2,tol)    

  # core eigenvalues 
  c<-eigen(KCD$C,symmetric=TRUE)$val 

  # objective function 
  lpSw<-function(w){    
    nu<-n*w/(1-w)+p1*p2+1 
    covKCD::lmvgamma((nu+n)/2,p1*p2)-covKCD::lmvgamma(nu/2,p1*p2) +
   .5*p1*p2*( nu*log(w) + n*log(1-w)) - 
   .5*(nu+n)*sum(log((1-w)*c+w))  
  }

  # maximize objective function 
  w<-optimize(lpSw,c(1e-6,1-1e-6),maximum=TRUE)$max 

  # Empirical Bayes posterior mean estimate
  SigmaHat<-(1-w)*S+w*KCD$K
  attr(SigmaHat,"w")<-w

  SigmaHat
}
