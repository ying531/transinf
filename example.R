
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
  return(c( weighted.orc$coefficients[2], sqrt(n)*sd(w.orc.influence)))
} 

weight.simu.fitted <- function(seed){
  set.seed(seed)
  data = gen.data() 
  dtrain = data.frame(cbind(data$X.org, data$Y.org)) 
  colnames(dtrain)[11] = "Y" 
  # fit the covariate shift 
  hat.ex.org = predict(regression_forest(data.frame(rbind(X.org, X.new)), c(rep(0,n), rep(1,m))), X.org)$predictions
  hat.wx.org = hat.ex.org * n / (m * (1-hat.ex.org))
  # plug in and weighted OLS
  weighted.lm = lm(Y~., data = dtrain, weights = hat.wx.org)
  w.influence = influence(weighted.lm)$coefficients[,2:3] 
  return(c( weighted.lm$coefficients[2], sqrt(n)*sd(w.influence)))
}
 

weight.simu <- function(seed){
  set.seed(seed)
  data = gen.data()  
  
  n = 1000
  m = 1000
  # our weighted procedure with correction 
  fold1 = sample(1:(n+m), floor((n+m)/2))  
  fold2 = setdiff(1:(n+m), fold1) 
  # cross-fitting the weights
  ws = rep(0, m+n) 
  df.Z = data.frame(rbind(data$X.org, data$X.new))
  label.Z = c(rep(0,n),rep(1,m)) 
  Z.rf.1 = regression_forest(df.Z[fold2,], label.Z[fold2])
  rf.pred.1 = predict(Z.rf.1, newdata = df.Z[fold1,])$predictions
  Z.rf.2 = regression_forest(df.Z[fold1,], label.Z[fold1])
  rf.pred.2 = predict(Z.rf.2, newdata = df.Z[fold2,])$predictions
  ws[fold1] = rf.pred.1 * n / ((1-rf.pred.1) * m)
  ws[fold2] = rf.pred.2 * n / ((1-rf.pred.2) * m)
  ws = ws[1:n]
  # fit weighted lm 
  lm.mdl = lm(Y~., data=dtrain, weights=ws)
  fitted.coef = coef(lm.mdl)[2]
  # infl.vals computes w_i * psi(D_i)
  infl.vals = n * influence(lm.mdl)$coefficients[,2] 
  # obtain estimate for psi(D_i)
  uw.infl.vals = infl.vals / ws   
  
  # conditional mean regression with grf
  fit.Zr = rep(0, n)
  fit.Zr.new = rep(0, m)  
  Zr.1 = regression_forest(data.frame(data$X.org[fold2[fold2<=n],]), uw.infl.vals[fold2[fold2<=n]], num.threads=1)
  Zr.2 = regression_forest(data.frame(data$X.org[fold1[fold1<=n],]), uw.infl.vals[fold1[fold1<=n]], num.threads=1)
  fit.Zr[fold1[fold1<=n]] = predict(Zr.1, data.frame(data$X.org[fold1[fold1<=n],]))$predictions
  fit.Zr[fold2[fold2<=n]] = predict(Zr.2, data.frame(data$X.org[fold2[fold2<=n],]))$predictions
  fit.Zr.new = (predict(Zr.1, data.frame(data$X.new))$predictions + predict(Zr.2, data.frame(data$X.new))$predictions)/2
  
  trans.coef = fitted.coef - mean(ws * fit.Zr, na.rm=TRUE) + mean(fit.Zr.new, na.rm=TRUE)
  
  trans.sd = sqrt(  sd(ws * (uw.infl.vals - fit.Zr), na.rm=TRUE)^2/n + 1/m * sd(fit.Zr.new, na.rm=TRUE)^2 )
  trans.cond.sd = sd(ws * (uw.infl.vals - fit.Zr), na.rm=TRUE)  
  
  return(c(trans.coef, trans.sd, trans.cond.sd))
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



