
gen.data <- function(){
  N = 10000
  n = 1000
  m = 1000
  p = 10
  X = matrix(runif(N*p), nrow=N, ncol=p)  
  noise = matrix(runif(N*2), nrow=N, ncol=2)  %*% chol(matrix(c(1,1/2,1/2,1), nrow=2))
  X[,3] = X[,1] + noise[,1]
  X[,4] = X[,1] + noise[,2]
  Y = X[,1] + abs(X[,1]-0.5) + X[,3] + rnorm(N) * 0.2 
  ex = (1 + pbeta(X[,1], shape1 = 3, shape2=2)) / 4
  pp = mean((1 + pbeta(runif(1000000), shape1 = 3, shape2=2)) / 4)
  wx = ex * (1-pp) / (pp * (1-ex)) 
  TT = rbinom(N, 1, ex)
  X.new = (X[TT==1,])[1:m,]
  Y.new = (Y[TT==1])[1:m]
  X.org = (X[TT==0,])[1:n,]
  Y.org = (Y[TT==0])[1:n]
  wx.org = (wx[TT==0])[1:n]
  wx.new = (wx[TT==1])[1:m]
  
  return(list("X.org"=X.org, "Y.org"=Y.org, "X.new"=X.new, "Y.new"=Y.new, "wx.org"=wx.org, "wx.new"=wx.new))
}

weight.simu.oracle <- function(seed){
  set.seed(seed)
  data = gen.data() 
  dtrain = data.frame(cbind(data$X.org, data$Y.org)) 
  colnames(dtrain)[11] = "Y" 
  weighted.orc = lm(Y~., data = dtrain, weights = data$wx.org)
  w.orc.influence = influence(weighted.orc)$coefficients[,2]
  return(c( weighted.orc$coefficients[2], sqrt(nrow(data$X.org))*sd(w.orc.influence)))
} 

weight.simu.fitted <- function(seed){
  set.seed(seed)
  data = gen.data() 
  dtrain = data.frame(cbind(data$X.org, data$Y.org)) 
  colnames(dtrain)[11] = "Y" 
  n = nrow(data$X.org)
  m = nrow(data$X.new)
  # fit the covariate shift and plug into weighted OLS
  hat.ex.org = predict(regression_forest(data.frame(rbind(data$X.org, data$X.new)), 
                                         c(rep(0,n), rep(1,m))), data$X.org)$predictions
  hat.wx.org = hat.ex.org * n / (m * (1-hat.ex.org)) 
  weighted.lm = lm(Y~., data = dtrain, weights = hat.wx.org) 
  return(c( weighted.lm$coefficients[2], sqrt(n)*sd(influence(weighted.lm)$coefficients[,2:3] )))
}
 

weight.simu.trans <- function(seed){
  set.seed(seed)
  data = gen.data()  
  dtrain = data.frame(cbind(data$X.org, data$Y.org)) 
  colnames(dtrain)[11] = "Y" 
  X.new = data.frame(data$X.new)
  lm.mdl = lm(Y~., data = dtrain)
  tlm.mdl = transfer(lm.mdl, df.new=X.new, alg='grf', verbose=FALSE)
  
  return(c(tlm.mdl$trans.coef[2], tlm.mdl$sup.w.std.err[2]/sqrt(nrow(data$X.org)), 
           tlm.mdl$trans.std.err[2]/sqrt(nrow(data$X.org))))
}


### the illustrating example 

# valid coverage using oracle weights
res.orc = sapply(1:1000, weight.simu.oracle)
mean( (1.0930 >= (res.orc[1,] -  res.orc[2,]*qnorm(0.975))) *
          + (1.0930 <= (res.orc[1,] +  res.orc[2,]*qnorm(0.975))) )

# undercover using fitted weights
res.fitted = sapply(1:1000, weight.simu.fitted)
mean( (1.0930 >= (res.fitted[1,] -  res.fitted[2,]*qnorm(0.975))) *
          + (1.0930 <= (res.fitted[1,] +  res.fitted[2,]*qnorm(0.975))) )

# valid coverage using our transductive estimator 
res.trans = sapply(1:1000, weight.simu.trans)
mean( (1.0930 >= (res.trans[1,] -  res.trans[2,]*qnorm(0.975))) *
          + (1.0930 <= (res.trans[1,] +  res.trans[2,]*qnorm(0.975))) )


# compare the lengths
mean.len.orc = mean(res.orc[2,])
mean.len.fitted = mean(res.fitted[2,])
mean.len.trans = mean(res.trans[2,])

c(mean.len.orc, mean.len.fitted, mean.len.trans)

# conditional transductive inference 
cat(mean(res.trans[3,]))



