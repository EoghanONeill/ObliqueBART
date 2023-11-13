
# ObliqueBART

<!-- badges: start -->
<!-- badges: end -->

The goal of ObliqueBART is to provide an implementation of Bayesian Additive Regression Trees with Oblique splits (i.e. hyperplane splits defined by linear combinations of covariates)

## Installation

``` r
library(devtools)
install_github("EoghanONeill/ObliqueBART")
```

## Example


``` r
library(ObliqueBART)
## example code


# example without hyperprior
library(dbarts)
library(dplyr)

# Simulate data -------------------------------------------
f <- function(x) {
  10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
    10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 200

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

# run BART --------------------------------------------------------
set.seed(99)
bartFit <- ObliqueBART(x, y, ntree = 10, nburn = 1000, npost = 5000)



cbind(apply(bartFit$y_hat,2,mean) , y)



ntest     <- 200

xtest  <- matrix(runif(ntest * 10), ntest, 10)
Eytest <- f(xtest)
ytest  <- rnorm(ntest, Eytest, sigma)


testres <- predict_mybart(bartFit, xtest,
               type = 'all')

cbind(apply(testres,2,mean) , ytest)



#### example of how the hyperprior can learn the important covariates


bartFit <- ObliqueBART(x, y, ntree = 10, nburn = 10000, npost = 10000,
                       coef_prior = "univariate_normal",
                       coef_hyperprior = "univariate_normal_betabinomial_theta_j",
                       coef_norm_hyperprior = "fixed")

testres <- predict_mybart(bartFit, xtest,
               type = 'all')

cbind(apply(testres,2,mean) , ytest)

bartFit$hyp_par_list$theta_vec



```

