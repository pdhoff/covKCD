covKCD
================
Peter Hoff
2022-07-27

### The Kronecker-core decomposition

Let
![Y\_1,\\ldots,Y\_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y_1%2C%5Cldots%2CY_n
"Y_1,\\ldots,Y_n") be a sample of ![p\_1\\times
p\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_1%5Ctimes%20p_2
"p_1\\times p_2") matrices, and let
![y\_1,\\ldots,y\_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_1%2C%5Cldots%2Cy_n
"y_1,\\ldots,y_n") be their vectorizations. The sample covariance matrix
![S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S
"S") of ![y\_1,\\ldots,
y\_n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_1%2C%5Cldots%2C%20y_n
"y_1,\\ldots, y_n") is a ![p\\times
p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%5Ctimes%20p
"p\\times p") positive semidefinite matrix, where
![p=p\_1p\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3Dp_1p_2
"p=p_1p_2"). The Kronecker-core decomposition
![(K,C)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28K%2CC%29
"(K,C)") of the covariance matrix
![S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S
"S") is a matrix factorization of the form   
![&#10;S = (K\_2\\otimes K\_1)^{1/2} C (K\_2\\otimes K\_1)^{1/2},
&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AS%20%3D%20%28K_2%5Cotimes%20K_1%29%5E%7B1%2F2%7D%20C%20%28K_2%5Cotimes%20K_1%29%5E%7B1%2F2%7D%2C%20%0A
"
S = (K_2\\otimes K_1)^{1/2} C (K_2\\otimes K_1)^{1/2}, 
")  
where

  - ![K\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_1
    "K_1") is a ![p\_1\\times
    p\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_1%5Ctimes%20p_1
    "p_1\\times p_1") matrix representing across-row covariance of the
    data matrices;  
  - ![K\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_2
    "K_2") is a ![p\_2\\times
    p\_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_2%5Ctimes%20p_2
    "p_2\\times p_2") matrix representing across-column covariance of
    the data matrices;
  - ![C](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C
    "C") is a core covariance matrix with zero average across-row and
    across-column covariance.

The covariance matrix ![K=K\_2\\otimes
K\_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K%3DK_2%5Cotimes%20K_1
"K=K_2\\otimes K_1") represents the “separable part” of
![S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S
"S"), whereas
![C](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C
"C") represents the “nonseparable part”. In particular, if
![S](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S
"S") is separable, than
![C=I\_p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%3DI_p
"C=I_p").

### Reference

Hoff, McCormack and Zhang (2022). [Core Shrinkage Covariance Estimation
for Matrix-variate Data.](https://arxiv.org/abs/2207.12484)

### Installation

``` r
# Development version 
devtools::install_github("pdhoff/covKCD")   

# CRAN-approved version 
install.packages("covKCD") 
```

### Basic usage

``` r
p1<-3 ; p2<-4 ; n<-10 

# Random sample of n independent standard normal p1xp2 matrices 
Y<-array(rnorm(n*p1*p2),dim=c(n,p1,p2)) 
dim(Y)
```

    ## [1] 10  3  4

``` r
# Sample covariance matrix 
S<-covKCD::mcov(Y) 
dim(S) 
```

    ## [1] 12 12

``` r
# KCD 
KCDS<-covKCD::covKCD(S,p1,p2) 

KCDS$K1 
```

    ##              [,1]       [,2]         [,3]
    ## [1,]  1.298126677 0.03912662 -0.002626438
    ## [2,]  0.039126618 0.94102020  0.107835507
    ## [3,] -0.002626438 0.10783551  0.795424559

``` r
KCDS$K2
```

    ##            [,1]        [,2]         [,3]         [,4]
    ## [1,]  0.8713638  0.25587722  0.102419374 -0.152234027
    ## [2,]  0.2558772  0.79230140  0.047605685 -0.284762704
    ## [3,]  0.1024194  0.04760568  0.844067547 -0.001949622
    ## [4,] -0.1522340 -0.28476270 -0.001949622  0.846966806

``` r
# Check decomposition 
K1h<-covKCD::msqrt(KCDS$K1) 
K2h<-covKCD::msqrt(KCDS$K2) 

range( kronecker(K2h,K1h) %*% KCDS$C %*% kronecker(K2h,K1h) - S )
```

    ## [1] -9.436896e-16  2.442491e-15

### Additional examples

  - [Simulation study from Hoff, McCormack and Zhang
    (2022).](https://www2.stat.duke.edu/~pdh10/Code/arXiv.2207.12484/)

  - [Speech recognition example from Hoff, McCormack and Zhang
    (2022).](https://www2.stat.duke.edu/~pdh10/Code/arXiv.2207.12484/)
