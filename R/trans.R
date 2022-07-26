#' Transductive inference for linear model 
#' 
#' Implement the transductive inference procedure for (weighted) linear regression, wrapping around lm() function 
#' @param formula The formula for lm() regression
#' @param data The dataframe for lm() regression 
#' @param param The coefficients to conduct transductive inference, can be a mixture of string names and integer indices
#' @param Z The dataframe for covariate shift attributes 
#' @param new.Z New data for the covariate shift attributes 
#' @param X The dataframe for conditioning attributes; need to be a subset of Z; if not given, we use Z
#' @param new.X New data for the conditioning attributes; need to be a subset of new.Z; if not given, we use new.Z
#' @param weights Optional, pre-specified covariate shift (weights); if not given, we automatically fit using grf package
#' @param alg Optional, a string for name of algorithm in fitting the conditional mean of influence functions, current options include `loess` and `grf`
#' @param random.seed Optional, random seed for sample splitting
#' @param other.params Optional, other parameters for the regression algorithm; can include span and degree for loess
#' @param folds Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting
#' @return Estimate for new conditional parameter, Estimate in default weighted OLS, standard errors and p-values for transductive inference, conditional confidence intervals for the new conditional parameter, standard errors and p-values for weighted OLS 
#' @examples 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1
#' Z = data.frame(X[,1:2])
#' new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
#' data = data.frame(cbind(X,Y))
#' lm.mdl = lm(Y~., data = data.frame(X))
#' trans.lm(Y~., data, c(1,"X1","X2"), Z, new.Z, alg="grf")
#' 
#' @export
trans.lm <- function(formula, data, param, Z, new.Z, 
                     X=NULL, new.X=NULL, weights=NULL, alg="loess",
                     random.seed=NULL,other.params=NULL,folds=NULL){
  n = nrow(data)
  m = nrow(new.Z)
  
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
  
  # sample splitting  
  if (is.null(folds)){
    fold1 = sample(1:(n+m), floor((n+m)/2)) 
  }else{
    fold1 = folds[[1]]
  } 
  fold2 = setdiff(1:(n+m), fold1) 
  
  ##=======================================##
  ##=== fit the weights if not provided ===## 
  ##=======================================## 
  if (is.null(weights)){
    ws = rep(0, m+n)
    
    df.Z = data.frame(rbind(as.matrix(Z), as.matrix(new.Z)))
    label.Z = c(rep(0, nrow(Z)), rep(1, nrow(new.Z)))
    # cross fitting 
    Z.rf.1 = regression_forest(df.Z[fold2,], label.Z[fold2])
    rf.pred.1 = predict(Z.rf.1, newdata = df.Z[fold1,])$predictions
    Z.rf.2 = regression_forest(df.Z[fold1,], label.Z[fold1])
    rf.pred.2 = predict(Z.rf.2, newdata = df.Z[fold2,])$predictions
    ws[fold1] = rf.pred.1 * n / ((1-rf.pred.1) * m)
    ws[fold2] = rf.pred.2 * n / ((1-rf.pred.2) * m)
    ws = ws[1:n]
  }else{
    ws = weights
  }
  
  ##=====================================##
  ##===== run weighted linear model =====## 
  ##=====================================##
  
  lm.mdl = lm(formula, data, weights=ws)
  
  names = param 
  for (i.par in 1:length(param)){
    if (!is.na(suppressWarnings(as.integer(param[i.par])))){
      names[i.par] = names(lm.mdl$coefficients)[as.integer(param[i.par])]
    }
  }
  
  fitted.coef = coef(lm.mdl)[names]
  # infl.vals computes w_i * psi(D_i)
  infl.vals = matrix(n * influence(lm.mdl)$coefficients[,names], ncol=length(param))
  # obtain estimate for psi(D_i)
  uw.infl.vals = diag(as.vector(1/ws)) %*% infl.vals  
  
  ##======================================##
  ##======= transductive inference =======## 
  ##======================================##
  
  trans.coefs = rep(0, length(param))
  trans.sds = rep(0, length(param))
  sup.sds = rep(0, length(param))
  
  for (i.par in 1:length(param)){ 
    ##=======================================##
    ##===== conditional mean regression =====## 
    fit.Zr = rep(0, n)
    fit.Zr.new = rep(0, m)
    # LOESS regression
    if (alg == 'loess'){
      loess.span = ifelse(is.null(other.params$span), 0.75, other.params$span)
      loess.deg = ifelse(is.null(other.params$degree), 2, other.params$degree) 
      # cross-fit nonparametric regression models
      Zr.1 = loess(uw.infl.vals[fold2[fold2<=n],i.par]~., data=data.frame(Z[fold2[fold2<=n],]),
                   span=loess.span, degree=loess.deg) 
      Zr.2 = loess(uw.infl.vals[fold1[fold1<=n],i.par]~., data=data.frame(Z[fold1[fold1<=n],]),
                   span=loess.span, degree=loess.deg) 
      fit.Zr[fold1[fold1<=n]] = predict(Zr.1, data.frame(Z[fold1[fold1<=n],]))
      fit.Zr[fold2[fold2<=n]] = predict(Zr.2, data.frame(Z[fold2[fold2<=n],]))
      fit.Zr.new = (predict(Zr.1, data.frame(new.Z)) + predict(Zr.2, data.frame(new.Z)))/2
    } 
    # random forest regression
    if (alg == 'grf'){
      # cross-fit nonparametric regression models
      Zr.1 = regression_forest(data.frame(Z[fold2[fold2<=n],]), 
                               uw.infl.vals[fold2[fold2<=n],i.par], num.threads=1)
      Zr.2 = regression_forest(data.frame(Z[fold1[fold1<=n],]), 
                               uw.infl.vals[fold1[fold1<=n],i.par], num.threads=1)
      fit.Zr[fold1[fold1<=n]] = predict(Zr.1, data.frame(Z[fold1[fold1<=n],]))$predictions
      fit.Zr[fold2[fold2<=n]] = predict(Zr.2, data.frame(Z[fold2[fold2<=n],]))$predictions
      fit.Zr.new = (predict(Zr.1, data.frame(new.Z))$predictions + predict(Zr.2, data.frame(new.Z))$predictions)/2
    }
    
    ##======================================##
    ##=== compute transductive estimator ===##  
    trans.coefs[i.par] = fitted.coef[i.par] - mean(ws * fit.Zr, na.rm=TRUE) + mean(fit.Zr.new, na.rm=TRUE)
    
    ##=======================================##
    ##== compute super-population variance ==##   
    sup.sds[i.par] = sqrt(sd(ws * (uw.infl.vals[,i.par] - fit.Zr), na.rm=TRUE)^2 + n/m * sd(fit.Zr.new, na.rm=TRUE)^2)
    
    ##===============================================##
    ##== compute conditional-transductive variance ==##   
    if (is.null(X) | ncol(X) == ncol(Z)){
      # If X=Z, then directly compute std
      trans.sds[i.par] = min(sup.sds[i.par], sd(ws * (uw.infl.vals[,i.par] - fit.Zr), na.rm=TRUE))
    }else{
      # cross-fit E[phi(D)|Z] on X in the new population
      fit.Xr = rep(0, m)
      # LOESS regression
      if (alg == 'loess'){
        loess.span = ifelse(is.null(other.params$span), 0.75, other.params$span)
        loess.deg = ifelse(is.null(other.params$degree), 2, other.params$degree) 
        # cross-fit nonparametric regression models
        Xr.1 = loess(fit.Zr.new[fold2[fold2>n]]~., data=data.frame(X.new[fold2[fold2>n],]),
                     span=loess.span, degree=loess.deg) 
        Xr.2 = loess(fit.Zr.new[fold1[fold1>n]]~., data=data.frame(X.new[fold1[fold1>n],]),
                     span=loess.span, degree=loess.deg) 
        fit.Xr[fold1[fold1>n]] = predict(Xr.1, data.frame(X.new[fold1[fold1>n],]))
        fit.Xr[fold2[fold2>n]] = predict(Xr.2, data.frame(X.new[fold2[fold2>n],])) 
      } 
      # random forest regression
      if (alg == 'grf'){
        # cross-fit nonparametric regression models
        Zr.1 = regression_forest(data.frame(X.new[fold2[fold2>n],]), 
                                 fit.Zr.new[fold2[fold2>n]], num.threads=1)
        Zr.2 = regression_forest(data=data.frame(X.new[fold1[fold1>n],]), 
                                 fit.Zr.new[fold1[fold1>n]], num.threads=1)
        fit.Xr[fold1[fold1>n]] = predict(Zr.1, data.frame(X.new[fold1[fold1>n],]))$predictions
        fit.Xr[fold2[fold2>n]] = predict(Zr.2, data.frame(X.new[fold2[fold2>n],]))$predictions 
      }
      # compute the variance
      trans.sds[i.par] = sqrt(min(sd(ws*uw.infl.vals[,i.par], na.rm=TRUE)^2, 
                             mean(ws^2*(uw.infl.vals[,i.par] - fit.Zr)^2, na.rm=TRUE)) + 
        n/m * mean((fit.Zr.new - fit.Xr)^2, na.rm=TRUE))
    }
    
    
  }
  
  trans.ci.lo = trans.coefs - qnorm(0.975) * trans.sds / sqrt(n)
  trans.ci.hi = trans.coefs + qnorm(0.975) * trans.sds / sqrt(n)
  trans.p_vals = 2 * pnorm(abs(trans.coefs) / (trans.sds/sqrt(n)), lower.tail = FALSE) 
  sup.p_vals = 2 * pnorm(abs(fitted.coef) / (sup.sds/sqrt(n)), lower.tail = FALSE)
  
  # print the results
  ret_table = cbind(trans.coefs, fitted.coef, trans.sds, trans.p_vals, sup.sds, sup.p_vals)
  colnames(ret_table) = c("Trans. Estimate", "W. OLS Estimate", 
                          "Trans. Std. Error", "Trans. Pr(>|z|)",
                          "Sup. W. Std. Error", "Sup. W. Pr(>|z|)")
  rownames(ret_table) = names
  
  cat("\n")
  cat("Summary of transductive inference in linear models")
  cat("\n\n")
  print(ret_table)
  cat("\n")
  
  
  invisible(list("trans.coef" = trans.coefs, "OLS.coef" = fitted.coef, 
                 "trans.std.err" = trans.sds, "trans.pval" = trans.p_vals,
                 "trans.ci.low" = trans.ci.lo, "trans.ci.upp" = trans.ci.hi,
                 "sup.w.std.err" = sup.sds, "sup.w.pval" = sup.p_vals))
   
}


#' Transductive inference for generalized linear model 
#' 
#' Implement the transductive inference procedure for (weighted) generalized linear regression, wrapping around glm() function 
#' @param formula The formula for glm() regression
#' @param family The family parameter for glm() regression
#' @param data The dataframe for glm() regression 
#' @param param The coefficients to conduct transductive inference, can be a mixture of string names and integer indices
#' @param Z The dataframe for covariate shift attributes 
#' @param new.Z New data for the covariate shift attributes 
#' @param X The dataframe for conditioning attributes; provide only when it differs from Z; need to be a subset of Z 
#' @param new.X New data for the conditioning attributes; provide only when it differs from new.Z; need to be a subset of new.Z 
#' @param weights Optional, pre-specified covariate shift (weights); if not given, we automatically fit using grf package
#' @param alg Optional, a string for name of algorithm in fitting the conditional mean of influence functions, current options include `loess` and `grf`
#' @param random.seed Optional, random seed for sample splitting
#' @param other.params Optional, other parameters for the regression algorithm; can include span and degree for loess
#' @param folds Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting
#' @return Estimate for new conditional parameter, Estimate in default weighted OLS, standard errors and p-values for transductive inference, conditional confidence intervals for the new conditional parameter, standard errors and p-values for weighted OLS 
#' @examples 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' logit.x = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + X[,1]**2 + rnorm(1000) * 0.1
#' Y = rbinom(n, 1, exp(logit.x)/(1+exp(logit.x)))
#' data = data.frame(cbind(X,Y))
#' Z = data.frame(X[,1:2])
#' new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1)
#' trans.glm(Y~., family='binomial', data, c(1,"V1","V2"), Z, new.Z, alg="grf")
#' 
#' @export
trans.glm <- function(formula, family, data, param, Z, new.Z, weights=NULL, alg="loess",
                     random.seed=NULL,other.params=NULL,folds=NULL){
  n = nrow(data)
  m = nrow(new.Z)
  
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
  
  # sample splitting  
  if (is.null(folds)){
    fold1 = sample(1:(n+m), floor((n+m)/2)) 
  }else{
    fold1 = folds[[1]]
  } 
  fold2 = setdiff(1:(n+m), fold1) 
  
  ##=======================================##
  ##=== fit the weights if not provided ===## 
  ##=======================================## 
  if (is.null(weights)){
    ws = rep(0, m+n)
    
    df.Z = data.frame(rbind(as.matrix(Z), as.matrix(new.Z)))
    label.Z = c(rep(0, nrow(Z)), rep(1, nrow(new.Z)))
    # cross fitting 
    Z.rf.1 = regression_forest(df.Z[fold2,], label.Z[fold2])
    rf.pred.1 = predict(Z.rf.1, newdata = df.Z[fold1,])$predictions
    Z.rf.2 = regression_forest(df.Z[fold1,], label.Z[fold1])
    rf.pred.2 = predict(Z.rf.2, newdata = df.Z[fold2,])$predictions
    ws[fold1] = rf.pred.1 * n / ((1-rf.pred.1) * m)
    ws[fold2] = rf.pred.2 * n / ((1-rf.pred.2) * m)
    ws = ws[1:n]
  }else{
    ws = weights
  }
  
  ##====================================##
  ##========= run weighted GLM =========## 
  ##====================================##
  
  glm.mdl = glm(formula, family = family, data, weights=ws)
  
  names = param 
  for (i.par in 1:length(param)){
    if (!is.na(suppressWarnings(as.integer(param[i.par])))){
      names[i.par] = names(glm.mdl$coefficients)[as.integer(param[i.par])]
    }
  }
  
  fitted.coef = coef(glm.mdl)[names]
  # infl.vals computes w_i * psi(D_i)
  infl.vals = matrix(n * influence(glm.mdl)$coefficients[,names], ncol=length(param))
  # obtain estimate for psi(D_i)
  uw.infl.vals = diag(as.vector(1/ws)) %*% infl.vals  
  
  ##======================================##
  ##======= transductive inference =======## 
  ##======================================##
  
  trans.coefs = rep(0, length(param))
  trans.sds = rep(0, length(param))
  sup.sds = rep(0, length(param))
  
  for (i.par in 1:length(param)){ 
    ##=======================================##
    ##===== conditional mean regression =====## 
    fit.Zr = rep(0, n)
    fit.Zr.new = rep(0, m)
    # LOESS regression
    if (alg == 'loess'){
      loess.span = ifelse(is.null(other.params$span), 0.75, other.params$span)
      loess.deg = ifelse(is.null(other.params$degree), 2, other.params$degree) 
      # cross-fit nonparametric regression models
      Zr.1 = loess(uw.infl.vals[fold2[fold2<=n],i.par]~., data=data.frame(Z[fold2[fold2<=n],]),
                   span=loess.span, degree=loess.deg) 
      Zr.2 = loess(uw.infl.vals[fold1[fold1<=n],i.par]~., data=data.frame(Z[fold1[fold1<=n],]),
                   span=loess.span, degree=loess.deg) 
      fit.Zr[fold1[fold1<=n]] = predict(Zr.1, data.frame(Z[fold1[fold1<=n],]))
      fit.Zr[fold2[fold2<=n]] = predict(Zr.2, data.frame(Z[fold2[fold2<=n],]))
      fit.Zr.new = (predict(Zr.1, data.frame(new.Z)) + predict(Zr.2, data.frame(new.Z)))/2
    }
    
    # random forest regression
    if (alg == 'grf'){
      # cross-fit nonparametric regression models
      Zr.1 = regression_forest(data.frame(Z[fold2[fold2<=n],]), uw.infl.vals[fold2[fold2<=n],i.par], num.threads=1)
      Zr.2 = regression_forest(data.frame(Z[fold1[fold1<=n],]), uw.infl.vals[fold1[fold1<=n],i.par], num.threads=1)
      fit.Zr[fold1[fold1<=n]] = predict(Zr.1, data.frame(Z[fold1[fold1<=n],]))$predictions
      fit.Zr[fold2[fold2<=n]] = predict(Zr.2, data.frame(Z[fold2[fold2<=n],]))$predictions
      fit.Zr.new = (predict(Zr.1, data.frame(new.Z))$predictions + predict(Zr.2, data.frame(new.Z))$predictions)/2
    }
    
    ##======================================##
    ##=== compute transductive estimator ===##  
    trans.coefs[i.par] = fitted.coef[i.par] - mean(ws * fit.Zr, na.rm=TRUE) + mean(fit.Zr.new, na.rm=TRUE)
    
    ##=======================================##
    ##== compute super-population variance ==##   
    sup.sds[i.par] = sqrt(sd(ws * (uw.infl.vals[,i.par] - fit.Zr), na.rm=TRUE)^2 + n/m * sd(fit.Zr.new, na.rm=TRUE)^2)
    
    ##===============================================##
    ##== compute conditional-transductive variance ==##   
    if (is.null(X) | ncol(X) == ncol(Z)){
      # If X=Z, then directly compute std
      trans.sds[i.par] = min(sup.sds[i.par], sd(ws * (uw.infl.vals[,i.par] - fit.Zr), na.rm=TRUE))
    }else{
      # cross-fit E[phi(D)|Z] on X in the new population
      fit.Xr = rep(0, m)
      # LOESS regression
      if (alg == 'loess'){
        loess.span = ifelse(is.null(other.params$span), 0.75, other.params$span)
        loess.deg = ifelse(is.null(other.params$degree), 2, other.params$degree) 
        # cross-fit nonparametric regression models
        Xr.1 = loess(fit.Zr.new[fold2[fold2>n]]~., data=data.frame(X.new[fold2[fold2>n],]),
                     span=loess.span, degree=loess.deg) 
        Xr.2 = loess(fit.Zr.new[fold1[fold1>n]]~., data=data.frame(X.new[fold1[fold1>n],]),
                     span=loess.span, degree=loess.deg) 
        fit.Xr[fold1[fold1>n]] = predict(Xr.1, data.frame(X.new[fold1[fold1>n],]))
        fit.Xr[fold2[fold2>n]] = predict(Xr.2, data.frame(X.new[fold2[fold2>n],])) 
      } 
      # random forest regression
      if (alg == 'grf'){
        # cross-fit nonparametric regression models
        Zr.1 = regression_forest(data.frame(X.new[fold2[fold2>n],]), 
                                 fit.Zr.new[fold2[fold2>n]], num.threads=1)
        Zr.2 = regression_forest(data=data.frame(X.new[fold1[fold1>n],]), 
                                 fit.Zr.new[fold1[fold1>n]], num.threads=1)
        fit.Xr[fold1[fold1>n]] = predict(Zr.1, data.frame(X.new[fold1[fold1>n],]))$predictions
        fit.Xr[fold2[fold2>n]] = predict(Zr.2, data.frame(X.new[fold2[fold2>n],]))$predictions 
      }
      # compute the variance
      trans.sds[i.par] = sqrt(min(sd(ws*uw.infl.vals[,i.par], na.rm=TRUE)^2, 
                                  mean(ws^2*(uw.infl.vals[,i.par] - fit.Zr)^2, na.rm=TRUE)) + 
                                n/m * mean((fit.Zr.new - fit.Xr)^2, na.rm=TRUE))
    }
    
  }
  
  trans.ci.lo = trans.coefs - qnorm(0.975) * trans.sds / sqrt(n)
  trans.ci.hi = trans.coefs + qnorm(0.975) * trans.sds / sqrt(n)
  trans.p_vals = 2 * pnorm(abs(trans.coefs) / (trans.sds/sqrt(n)), lower.tail = FALSE) 
  sup.p_vals = 2 * pnorm(abs(fitted.coef) / (sup.sds/sqrt(n)), lower.tail = FALSE)
  
  # print the results
  ret_table = cbind(trans.coefs, fitted.coef, trans.sds, trans.p_vals, sup.sds, sup.p_vals)
  colnames(ret_table) = c("Trans. Estimate", "W. OLS Estimate", 
                          "Trans. Std. Error", "Trans. Pr(>|z|)",
                          "Sup. W. Std. Error", "Sup. W. Pr(>|z|)")
  rownames(ret_table) = names
  
  cat("\n")
  cat("Summary of transductive inference in generalized linear models")
  cat("\n\n")
  print(ret_table)
  cat("\n")
  
  
  invisible(list("trans.coef" = trans.coefs, "OLS.coef" = fitted.coef, 
                 "trans.std.err" = trans.sds, "trans.pval" = trans.p_vals,
                 "trans.ci.low" = trans.ci.lo, "trans.ci.upp" = trans.ci.hi,
                 "sup.w.std.err" = sup.sds, "sup.w.pval" = sup.p_vals))
  
}

