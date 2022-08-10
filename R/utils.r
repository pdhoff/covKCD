#' Covariance matrix to covariance array
#' 
#' Reshape a covariance matrix to a covariance array. 
#' 
#' @param S a covariance matrix of dimension (p1p2)*(p1p2).
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' 
#' @return a four-way array where entry i1,j1,i2,j2 gives the covariance 
#' between element i1,j1 and element i2,j2 of a random matrix.
#' 
#' @author Peter Hoff
#'
#' @examples
#' p1<-4 ; p2<-7 ; p<-p1*p2 
#'  
#' S<-rWishart(1,p,diag(p))[,,1] 
#' A<-cm2ca(S,p1,p2) 
#' range(S-ca2cm(A)) 
#' 
#' @export
cm2ca<-function(S,p1,p2){
  A<-array(dim=c(p1,p2,p1,p2)) 
  for(j1 in 1:p2){ for(j2 in 1:p2){
    A[,j1,,j2]<-S[p1*(j1-1)+1:p1,p1*(j2-1)+1:p1]
  }} 
  A
}


#' Covariance array to covariance matrix
#' 
#' Reshape a covariance array to a covariance matrix.
#' 
#' @param A a covariance array of dimension p1*p2*p1*p2. 
#'  
#' @return a p1*p2 by p1*p2 covariance matrix.
#'
#' @author Peter Hoff
#'
#' @examples
#' p1<-4 ; p2<-7 ; p<-p1*p2 
#' 
#' S<-rWishart(1,p,diag(p))[,,1] 
#' A<-cm2ca(S,p1,p2) 
#' range(S-ca2cm(A)) 
#' 
#' @export
ca2cm<-function(A){ 
  p1<-dim(A)[1] ; p2<-dim(A)[2] 
  S<-matrix(0,p1*p2,p1*p2) 
  for(j1 in 1:p2){ for(j2 in 1:p2){
    S[p1*(j1-1)+1:p1,p1*(j2-1)+1:p1]<-A[,j1,,j2]
  }} 
  S
}


#' Matrix-variate covariance matrix
#' 
#' Compute the covariance matrix of a sample of data matrices.
#' 
#' @param Y a numeric n*p1*p2 data array corresponding to n data matrices 
#' of dimension p1*p2.
#' @param use a character string giving method for dealing with missing 
#' values, fed to the \code{\link[stats]{cov}} function.  
#'
#' @return a p1*p2 by p1*p2 sample covariance matrix of the n vectorized data
#' matrices.
#'  
#' @author Peter Hoff 
#' 
#' @examples 
#' p1<-4 ; p2<-3 ; n<-200
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
#' image(S)     
#'
#' plot(S,kronecker(Sig2,Sig1)) ; abline(0,1) 
#'
#' @import stats
#' @export
mcov<-function(Y,use="everything"){ cov(t(apply(Y,1,c)),use=use) }


#' Symmetric square root of a matrix
#' 
#' Compute the symmetric square root of a matrix.
#' 
#' @param M a positive semidefinite matrix. 
#'
#' @return a positive semidefinite matrix. 
#' 
#' @author Peter Hoff 
#' 
#' @examples
#' S<-rWishart(1,5,diag(5))[,,1]
#' S 
#' Sh<-msqrt(S)
#' Sh%*%Sh 
#' 
#' @export
msqrt<-function(M)
{
  tmp<-eigen(M,symmetric=TRUE)
  tmp$vec%*%(t(tmp$vec)*sqrt(tmp$val))
} 

#' Inverse symmetric square root of a matrix
#' 
#' Compute the inverse of the symmetric square root of a matrix.
#' 
#' @param M a positive definite matrix.  
#' 
#' @return a positive definite matrix.
#'
#' @author Peter Hoff 
#' 
#' @examples
#' S<-rWishart(1,5,diag(5))[,,1]
#' solve(S) 
#' iSh<-msqrtInv(S)
#' iSh%*%iSh 
#' 
#' @export
msqrtInv<-function(M)
{
  tmp<-eigen(M,symmetric=TRUE)
  tmp$vec%*%(t(tmp$vec)*sqrt(1/tmp$val))
} 


