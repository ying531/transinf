# transinf
R package that implements conditional transductive inference of parameters in lm and glm models. Details of conditional and transductive inference is in [this paper](https://arxiv.org/abs/2104.04565). The current implementation uses the influence function from `lm()` and `glm()` in R and two-fold cross-fitting for conditional standard error estimation. 

## Installation

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("ying531/transinf")`.

## Why using our method

We illustrate the necessity of using our method when there is a covariate shift and covariates from the target distribution are observed. Before talking about conditional inference, we look at the conventional setting where the target is a super-population parameter. We will use the example of ordinary least squares (OLS) under covariate shifts. 



In our simulation setting, there is an unknown covariate shift on all the regression variables from the training to the target distribution. The covariate shift gives 'weights' to the training data. We have access to the original dataset of size `n=1000`, and the covariates of the new dataset of size `m=1000`. From Monte Carlo, the true super-population parameter under the target distribution is `1.0930`.



All the codes are in `example.R`, where `gen.data()` generates synthetic datasets, `weight.simu.oracle(seed)` runs weighted OLS with true weights, `weight.simu.fitted(seed)` runs weighted OLS with fitted weights, `weight.simu(seed)` runs our procedure. To run all the experiments below, simply load the file `example.R`. 



If the weights are known, one could run weighted OLS and get valid coverage, as follows: 

```R
> res.orc = sapply(1:1000, weight.simu.oracle)
> mean( (1.0930 >= (res.orc[1,] -  res.orc[2,]*qnorm(0.975))) *
+         (1.0930 <= (res.orc[1,] +  res.orc[2,]*qnorm(0.975))) )
[1] 0.953
```



In practice the weights are typically unknown. The common way is to run weighted OLS with fitted weights, i.e., first fit a weight (perhaps using machine learning), and then run OLS with weights on the training set. The following codes uses `regression_forest` from the `grf` package to fit the weights and then construct confidence intervals with weighted OLS. 

```R
> res.fitted = sapply(1:1000, weight.simu.fitted)
> mean( (1.0930 >= (res.fitted[1,] -  res.fitted[2,]*qnorm(0.975))) *
+         (1.0930 <= (res.fitted[1,] +  res.fitted[2,]*qnorm(0.975))) )

```







## Helper

##### Usage


```
trans.lm(
  trans.lm(
  formula,
  data,
  param,
  Z,
  new.Z,
  X = NULL,
  new.X = NULL,
  weights = NULL,
  alg = "loess",
  random.seed = NULL,
  other.params = NULL,
  folds = NULL
)
```

This function builds upon the `lm()` or `glm()` functions, and have similar usage as them. It prints the summary of transductive inference on the specified coefficients, and returns a list of results as described in the following. 

| `formula`      | The formula for lm() regression                              |
| -------------- | ------------------------------------------------------------ |
| `data`         | The dataframe for lm() regression                            |
| `param`        | The coefficients to conduct transductive inference, can be a mixture of string names and integer indices |
| `Z`            | The dataframe of the covariate shift attributes              |
| `new.Z`        | New data for the covariate shift attributes                  |
| `X`            | The dataframe of the conditioning attributes; provide only when it differs from Z; need to be a subset of Z |
| `new.X`        | New data for the conditioning attributes; provide only when it differs from new.Z; need to be a subset of new.Z |
| `weights`      | Optional, pre-specified covariate shift (weights); if not given, we automatically fit using grf package |
| `alg`          | Optional, a string for name of algorithm in fitting the conditional mean of influence functions, current options include 'loess' and 'grf' |
| `random.seed`  | Optional, random seed for sample splitting                   |
| `other.params` | Optional, other parameters for the regression algorithm; can include span and degree for loess |
| `folds`        | Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting |

If `weights` is not given, or you would like to use `alg = 'grf'` as the regressor, the R package `grf` is required to be installed.

| Output          | Description                                                  |
| --------------- | ------------------------------------------------------------ |
| `trans.coef`    | Estimate for new conditional parameter                       |
| `OLS.coef`      | Estimate in default weighted OLS                             |
| `trans.std.err` | Standard errors for transductive inference                   |
| `trans.pval`    | P-values for transductive inference (testing whether the new conditional parameter is zero) |
| `trans.ci.low`  | Lower 0.95-confidence bound for the new (transductive) conditional parameter |
| `trans.ci.upp`  | Upper 0.95-confidence bound for the new (transductive) conditional parameter |
| `sup.w.std.err` | Standard errors for super-population inference using our method |
| `sup.w.pval`    | P-values for super-population parameter (testing for whether the new super-population parameter is zero, built around `trans.coef`) |


## Examples

### Transductive inference for linear models

The following example works out transductive inference of linear regression coefficients (setting `param=1` selects the intercept) for a well-specified linear model. We use `grf` to fit the conditional mean functions of the influence function.  By default, the covariate shifts hold for Z, and we conduct Z-conditional inference.



```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1) 
> lm.mdl = lm(Y~., data = data.frame(X))
> trans.lm(Y~., data = data.frame(X), c(1,"X1","X2"), Z, new.Z, X=NULL, new.X=NULL, alg="grf")

Summary of transductive inference in linear models

            Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error Sup. W. Pr(>|z|)
(Intercept)     0.001918684         0.1584812       0.7018331          0.1642096        0.7117612
X1              0.990066950         0.1557365       0.0000000          0.1600513        0.0000000
X2              1.998203754         0.1542064       0.0000000          0.1578653        0.0000000

```

In the above summary, `Trans. Estimate` is the estimator for the new conditional parameter, `Trans. Std. Err` is the estimated standard error for inferring the new conditional parameter, and `Trans. Pr(>|z|)` is the p-value for testing whether the new conditional parameter is zero. `Sup. W. Std. Error` and `Sup. W. Pr(>|z|)` are the standard error and p-value for the super-population parameter built around `Trans. Estimate`. We note that this is different from weighted OLS. 



The following example conducts conditional inference for a misspecified linear model. The regression algorithm is `grf` and we focus on the coefficients for the first coefficient (the intercept), `"X1"`, and `"X2"`.

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1) 
> lm.mdl = lm(Y~., data = data.frame(X))
> trans.lm(Y~., data = data.frame(X), c(1,"X1","X2"), Z, new.Z, alg="grf")

Summary of transductive inference in linear models

            Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error Sup. W. Pr(>|z|)
(Intercept)       0.6763192          2.105422    3.049578e-24           2.215640     4.783930e-22
X1                0.9844177          4.776830    7.178563e-11           4.786958     7.868624e-11
X2                2.1056581          3.048357   8.963295e-106           3.167178     3.956157e-98
```



#### Transductive inference for generalized linear models 

The following example works out transductive inference for parameters from a well-specified logistic model, and the weights are estimated. The regression algorithm is `grf`. We focus on the coefficients for the second variable (`"X1"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x))) 
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> trans.glm(Y~., family='binomial', data = data.frame(X), "X1", Z, new.Z, alg="grf")

Summary of transductive inference in generalized linear models

   Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error Sup. W. Pr(>|z|)
X1        1.136582          7.211347    6.226018e-07           7.606262     2.297764e-06
```

This is an example of conditional inference for parameters from a misspecified logistic regression model. The regression algorithm is `grf`. We focus on the coefficients for variable `"X3"` and the third variable (`"X2"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x))) 
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> trans.glm(Y~., family='binomial', data = data.frame(X), c("X3",3), Z, new.Z, alg="grf")

Summary of transductive inference in generalized linear models

   Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error Sup. W. Pr(>|z|)
X3        2.633818          7.511542    1.433210e-28           7.901520     5.599850e-26
X2        1.474090          6.071821    1.625554e-14           6.296418     1.327544e-13
```