# transinf
R package that implements conditional transductive inference of parameters in lm and glm models. Details of conditional inference is in [this paper](https://arxiv.org/abs/2104.04565). The current implementation uses the influence function from `lm()` and `glm()` in R and two-fold cross-fitting for conditional standard error estimation. 

## Installation

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("ying531/transinf")`.

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
| `Z`            | The dataframe of the conditioning set                        |
| `new.Z`        | New data for the conditioning variables                      |
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
| `sup.w.std.err` | Standard errors for weighted OLS                             |
| `sup.w.pval`    | P-values for weighted OLS (testing for whether the new super-population parameter is zero, not necessarily valid with estimated weights) |


## Examples

### Transductive inference for linear models

The following example works out transductive inference of linear regression coefficients (setting `param=1` selects the intercept) for a well-specified linear model. We use `grf` to fit the conditional mean functions of the influence function. 



```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> data = data.frame(cbind(X,Y))
> lm.mdl = lm(Y~., data = data.frame(X))
> trans.lm(Y~., data, c(1,"X1","X2"), Z, new.Z, alg="grf")

Summary of transductive inference in linear models

            Trans. Estimate W. OLS Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error
(Intercept)   -8.716312e-18   -8.695553e-18      1.030577e-17    1.38906e-157       1.030577e-17
X1             2.633077e-16    2.570466e-16      1.493904e-16     0.00000e+00       1.493904e-16
X2             4.457330e-16    4.385998e-16      3.015352e-16     0.00000e+00       3.015352e-16
            Sup. W. Pr(>|z|)
(Intercept)    7.633658e-157
X1              0.000000e+00
X2              0.000000e+00

```

In the above summary, `Trans. Estimate` is the estimator for the new conditional parameter, `W. OLS Estimate` is the weighted OLS estimator, `Trans. Std. Err` is the estimated standard error for inferring the new conditional parameter, and `Trans. Pr(>|z|)` is the p-value for testing whether the new conditional parameter is zero. `Sup. W. Std. Error` and `Sup. W. Pr(>|z|)` are the standard error and p-value for standard super-population inference from weighted OLS. 



The following example conducts conditional inference for a misspecified linear model. In this case, the inference for the conditional parameter can be different from that for the super-population parameter (see the `cond.std.err` for conditional inference and `std.err` for super-population inference). The regression algorithm is `grf` and we focus on the coefficients for the first coefficient, `"X1"`, and `"X2"`.

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> data = data.frame(cbind(X,Y))
> lm.mdl = lm(Y~., data = data.frame(X))
> trans.lm(Y~., data, c(1,"X1","X2"), Z, new.Z, alg="grf")

Summary of transductive inference in linear models

            Trans. Estimate W. OLS Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error
(Intercept)    3.298491e-18               0      1.221473e-16       0.3931339       1.221473e-16
X1            -3.074325e-18               0      1.803963e-16       0.5899439       1.803963e-16
X2             1.438591e-18               0      3.512461e-16       0.8969488       3.512461e-16
            Sup. W. Pr(>|z|)
(Intercept)                1
X1                         1
X2                         1
```



#### Transductive inference for generalized linear models 

The following example works out transductive inference for parameters from a well-specified logistic model, and the weights are estimated. The regression algorithm is `grf`. We focus on the coefficients for the second variable (`"V1"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
> data = data.frame(cbind(X,Y))
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> trans.glm(Y~., family='binomial', data, "V1", Z, new.Z, alg="grf")

Summary of transductive inference in generalized linear models

   Trans. Estimate W. OLS Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error Sup. W. Pr(>|z|)
V1       0.8393124       0.9987003          5.574298    1.922674e-06           5.574298     1.465214e-08
```

This is an example of conditional inference for parameters from a misspecified logistic regression model. The regression algorithm is `loess`. We focus on the coefficients for variable `"V3"` and the third variable (`"V2"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
> data = data.frame(cbind(X,Y))
> Z = data.frame(X[,1:2])
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> trans.glm(Y~., family='binomial', data, c("V3",3), Z, new.Z, alg="grf")

Summary of transductive inference in generalized linear models

   Trans. Estimate W. OLS Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. W. Std. Error Sup. W. Pr(>|z|)
V3        3.045811        2.809102          6.767647    5.801914e-46           6.767647     2.339215e-39
V2        1.589582        1.612555          5.259940    1.217242e-21           5.259940     3.176103e-22
```