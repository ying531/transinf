# transinf
R package that implements conditional transductive inference for lm and glm models. Details of conditional and transductive inference is in [this paper](https://arxiv.org/abs/2104.04565). The current implementation uses the influence function from `lm()` and `glm()` in R and two-fold cross-fitting for conditional standard error estimation. 

## Installation

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("ying531/transinf")`.

## Basic usage

Suppose one would like to transfer the linear model Y~X to a new population, and there is a covariate shift on attributes Z, a subset of the original covariates. Let `Z.new` store the new attributes, `data` be the original dataset containing Y and X (*column names for the attributes in Z should be consistent in `Z.new` and `data`*). Then simply run the original linear model with `lm()`. The `transfer` function in our package takes the `lm()` object as input and transfer it to the new population. 

```R
lm.mdl = lm(Y~., data = data) 
tlm.mdl = transfer(lm.mdl, newdata = Z.new)
```

By default, the `transfer()` function transfers all the coefficients in the `lm()` model. The usage is the same for `glm()` models. Below is an example for logistic regression. 

```R
glm.mdl = lm(Y~., family = 'binomial', data = data) 
tglm.mdl = transfer(glm.mdl, newdata = Z.new)
```

For details on the statistical inference procedure and theory, see [our paper](https://arxiv.org/abs/2104.04565). For details on the usage, see the Documentation and Examples parts below. 

## Why our method

We illustrate the necessity of using our method when there is a covariate shift and covariates from the target distribution are observed. Before talking about conditional inference, we look at the conventional setting where the target is a super-population parameter. We will use the example of ordinary least squares (OLS) under covariate shifts. 



In our simulation setting, there is an unknown covariate shift on all the regression variables from the training to the target distribution. The covariate shift gives 'weights' to the training data. We have access to the original dataset of size `n=1000`, and the covariates of the new dataset of size `m=1000`. From Monte Carlo, the true super-population parameter under the target distribution is `1.0930`.



All the codes are in `example.R`, where `gen.data()` generates synthetic datasets, `weight.simu.oracle(seed)` runs weighted OLS with true weights, `weight.simu.fitted(seed)` runs weighted OLS with fitted weights, `weight.simu(seed)` runs our procedure (an equivalent implementation as the package). To run all the experiments below, simply load the file `example.R`. 



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
[1] 0.870
```

There is severe undercover. This is because the weights often cannot be estimated accurately, incurring a bias which is much larger than the asymptotic variance. 



Our method builds confidence intervals around an estimator with a correction term. This allows us to conduct root-n inference when the weights are estimated at a lower rate. The weights are still estimated with `regression_forest` from the `grf` package. 

```R
> res.trans = sapply(1:1000, weight.simu.trans)
> mean( (1.0930 >= (res.trans[1,] -  res.trans[2,]*qnorm(0.975))) *
+         (1.0930 <= (res.trans[1,] +  res.trans[2,]*qnorm(0.975))) )
[1] 0.946
```

We indeed obtain valid coverage with moderate sample sizes. Indeed, our confidence intervals are even shorter than the previous two methods (without correction). We note that this is not a general rule: our confidence intervals might be longer if the size of the new datset, `m`, is smaller.

```R
mean.len.orc = mean(res.orc[2,])
mean.len.fitted = mean(res.fitted[2,])
mean.len.trans = mean(res.trans[2,])

c(mean.len.orc, mean.len.fitted, mean.len.trans)

[1] 0.04468186 0.04383533 0.04223867
```



Furthermore, if you would like to infer parameters about a finite population with observed covariates, turning to conditional parameter allows for shorter confidence intervals and conditional validity. 

```R
cat(mean(res.trans[3,]))

[1] 0.03998838
```





## Documentation


```
transfer( 
  object,
  newdata, 
  cond.newdata = NULL,
  param = NULL,
  wts = NULL,
  alg = "loess",
  random.seed = NULL,
  other.params = NULL,
  folds = NULL, 
  verbose = TRUE
)
```

This function takes a fitted `lm()` or `glm()` object as input, and transfers it to a new population.  You also need to specify `df.new`, the covariate shift attributes in the new population. **Please make sure the column names are consistent across the data used in lm() or glm() and in `df.new`**.

| Arguments      | Description                                                  |
| -------------- | ------------------------------------------------------------ |
| `object`    | An lm() or glm() object that fits the original regression    |
| `newdata` | Dataframe for covariate shift attributes for the new population |
| `cond.newdata` | Dataframe for the new conditioning set; default to be `df.new` if not provided; can be a subset of `df.new` |
| `param`        | The coefficients to conduct transductive inference; default to be all the original coefficients in mdl if not provided; can be a mixture of string names and integer indices |
| `wts`      | Optional, pre-specified covariate shift (weights); if not given, we automatically fit using grf package |
| `alg`          | Optional, a string for name of algorithm in fitting the conditional mean of influence functions, current options include 'loess' and 'grf' |
| `random.seed`  | Optional, random seed for sample splitting                   |
| `other.params` | Optional, other parameters for the regression algorithm; can include span and degree for loess |
| `folds`        | Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting |
| `verbose`        | Optimal, TRUE by default; print the summary if verbose==TRUE |

If `weights` is not given, or you would like to use `alg = 'grf'` as the regressor, the R package `grf` is required to be installed.

| Output          | Description                                                  |
| --------------- | ------------------------------------------------------------ |
| `trans.coef`    | Estimate for new conditional parameters                      |
| `trans.std.err` | Standard errors for transductive inference                   |
| `trans.pval`    | P-values for transductive inference (testing whether the new conditional parameter is zero) |
| `trans.ci.low`  | Lower 0.95-confidence bound for the new (transductive) conditional parameter |
| `trans.ci.upp`  | Upper 0.95-confidence bound for the new (transductive) conditional parameter |
| `sup.std.err` | Standard errors for super-population inference using our method |
| `sup.pval`    | P-values for super-population parameter (testing for whether the new super-population parameter is zero) using our method, built around `trans.coef`) |
| `summary`    | Summary table of the model fitting results; the printed result for `verbose=TRUE` |



## Examples

#### Transductive inference for linear models

The following example works out transductive inference of linear regression coefficients (setting `param=1` selects the intercept) for a well-specified linear model. We use `grf` to fit the conditional mean functions of the influence function.  By default, the covariate shifts hold for the attributes in Z (containing 'X1' and 'x2'), and we conduct Z-conditional inference. 



```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1 
> dat = data.frame(cbind(Y, data.frame(X)))
> colnames(dat)[1] = "Y" 
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1) 
> colnames(new.Z) = c("X1", "X2")
> lm.mdl = lm(Y~., data = dat)
> transfer(lm.mdl, newdata=new.Z, param=c(1,"X1","X2"), alg='grf') 

Summary of transductive inference in linear models

            Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
(Intercept)    -0.005514753         0.1458162       0.2317085       0.1461592     0.2328055
X1              0.994471361         0.2164838       0.0000000       0.2172301     0.0000000
X2              1.994978628         0.2250008       0.0000000       0.2254147     0.0000000

```

In the above summary, `Trans. Estimate` is the transductive estimator for the new conditional parameter, `Trans. Std. Err` is the estimated standard error for inferring the new conditional parameter, and `Trans. Pr(>|z|)` is the p-value for testing whether the new conditional parameter is zero. `Sup. Std. Error` and `Sup. Pr(>|z|)` are the standard error and p-value for the super-population parameter built around `Trans. Estimate`. We note that this is different from weighted OLS. 



The following example conducts conditional inference for a misspecified linear model. The regression algorithm is `grf` and we focus on the coefficients for the first coefficient (the intercept), `"X1"`, and `"X2"`.

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> dat = data.frame(cbind(Y, data.frame(X)))
> colnames(dat)[1] = "Y" 
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1) 
> colnames(new.Z) = c("X1", "X2")
> lm.mdl = lm(Y~., data = data.frame(X))
> transfer(lm.mdl, newdata=new.Z, param=c(1,"X1","X2"), alg="grf")

Summary of transductive inference in linear models

            Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
(Intercept)       0.3588687         0.2018726               0       0.4776382 8.773232e-125
X1                1.0103193         0.5844082               0       1.0553311 2.520910e-201
X2                2.0365865         0.4261175               0       0.7712821  0.000000e+00
```



#### Transductive inference for generalized linear models 

The following example works out transductive inference for parameters from a well-specified logistic model, and the weights are estimated. The regression algorithm is `loess`, with default `span=0.75` and `degree=2`. We focus on the coefficients for the second variable (`"X1"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
> Y = rbinom(1000, 1, exp(logit.x)/(1+exp(logit.x))) 
> dat = data.frame(cbind(Y, data.frame(X)))
> colnames(dat)[1] = "Y"   
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
> glm.mdl = glm(Y~., family='binomial', data=dat)
> transfer(glm.mdl, newdata=new.Z, param=c("X1"), alg='loess')

Summary of transductive inference in generalized linear models

   Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
X1        1.041567          7.403832    8.640142e-06        7.415148  8.917361e-06
```

This is an example of conditional inference for parameters from a misspecified logistic regression model. The regression algorithm is `grf`. We focus on the coefficients for variable `"X3"` and the third variable (`"X2"`).

```R
> X = matrix(rnorm(1000*10), nrow=1000)
> logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
> Y = rbinom(1000, 1, exp(logit.x)/(1+exp(logit.x))) 
> dat = data.frame(cbind(Y, data.frame(X)))
> colnames(dat)[1] = "Y"  
> new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1) 
> glm.mdl = glm(Y~., family='binomial', data=dat)
> transfer(glm.mdl, newdata=new.Z, param=c("X3",3), alg='grf')

Summary of transductive inference in generalized linear models

   Trans. Estimate Trans. Std. Error Trans. Pr(>|z|) Sup. Std. Error Sup. Pr(>|z|)
X3        3.150276          7.151847    4.202572e-44        7.202032  1.627970e-43
X2        2.021570          7.933163    7.737493e-16        7.981118  1.148373e-15

```

