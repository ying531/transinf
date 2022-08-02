#' Transductive inference for (generalized) linear models
#' 
#' Implement the transductive inference procedure for (weighted) linear and generalized linear regression, wrapping around an lm() or glm() object
#' @param object An lm() or glm() object that fits the original regression 
#' @param newdata Dataframe for covariate shift attributes for the new population
#' @param cond.newdata Optional, dataframe for the new conditioning set; default to be newdata if not provided; can be a subset of newdata
#' @param param The coefficients to conduct transductive inference; default to be all the original coefficients in object if not provided; can be a mixture of string names and integer indices 
#' @param wts Optional, pre-specified covariate shift (weights); if not given, we automatically fit using grf package
#' @param alg Optional, a string for name of algorithm in fitting the conditional mean of influence functions, current options include `loess` and `grf`
#' @param random.seed Optional, random seed for sample splitting
#' @param other.params Optional, a list of other parameters for the regression algorithm; can include span and degree for loess
#' @param folds Optional, a list of two folds of indices for sample splitting; can be useful to control sample splitting
#' @param verbose Optimal, TRUE by default; print the summary if verbose==TRUE
#' @return Point transductive estimator, standard errors and p-values for transductive inference, conditional confidence intervals for the new conditional parameters, standard errors and p-values for super-population parameter
#' @examples 
#' X = matrix(rnorm(1000*10), nrow=1000)
#' Y = X %*% matrix(c(1,2,3,rep(0,10-3)), ncol=1) + rnorm(1000) * 0.1 
#' dat = data.frame(cbind(Y, data.frame(X)))
#' colnames(dat)[1] = "Y"
#' lm.mdl = lm(Y~., data = dat)
#' new.Z = data.frame(matrix(runif(500*2), nrow=500)*2-1) 
#' colnames(new.Z) = c("X1", "X2")
#' transfer(lm.mdl, newdata=new.Z)
#' transfer(lm.mdl, newdata=new.Z, param=c(1,"X1","X2")) 
#' 
#' @export
transfer <- function(object, newdata, cond.newdata=NULL, 
                     param=NULL,wts=NULL,alg="loess",
                     random.seed=NULL,other.params=NULL,folds=NULL,verbose=TRUE){
  org.data = eval(getCall(object)$data)#object$model
  org.formula = eval(object$call[[2]])
  n = nrow(org.data)
  m = nrow(newdata)
  
  if (!is.null(random.seed)){
    set.seed(random.seed)
  }
   
  if (is.null(param)){
    param = 1:length(object$coefficients)
  }
  
  # sample splitting  
  if (is.null(folds)){ fold1 = sample(1:(n+m), floor((n+m)/2)) 
  }else{ fold1 = folds[[1]] } 
  fold2 = setdiff(1:(n+m), fold1) 
  
  # extract original Z 
  Z = org.data[,colnames(newdata)] 
  
  
  ##=======================================##
  ##=== fit the weights if not provided ===## 
  ##=======================================## 
  ws = rep(0, m+n)
  if (is.null(wts)){ 
    df.Z = data.frame(rbind(as.matrix(Z), as.matrix(newdata)))
    label.Z = c(rep(0, nrow(Z)), rep(1, nrow(newdata)))
    # cross fitting 
    Z.rf.1 = regression_forest(df.Z[fold2,], label.Z[fold2])
    rf.pred.1 = predict(Z.rf.1, newdata = df.Z[fold1,])$predictions
    Z.rf.2 = regression_forest(df.Z[fold1,], label.Z[fold1])
    rf.pred.2 = predict(Z.rf.2, newdata = df.Z[fold2,])$predictions
    ws[fold1] = rf.pred.1 * n / ((1-rf.pred.1) * m)
    ws[fold2] = rf.pred.2 * n / ((1-rf.pred.2) * m)
    ws = ws[1:n]
  }else{
    ws = wts
  } 
  
  ##=====================================##
  ##===== run weighted (generalized) linear model =====## 
  ##=====================================## 
  if (class(object)[1]=='lm'){
    wlm.mdl = lm.internal(formula=org.formula, data=org.data, weights=ws)
  }else if (class(object)[1]=='glm'){
    wlm.mdl = glm.internal(formula=org.formula, data=org.data, 
                           family=object$family$family, weights=ws)
  }else{
    cat("Not lm() or glm() model!")
    cat("\n")
    return()
  }
  
  # turn indices to coefficient names
  names = param 
  for (i.par in 1:length(param)){
    if (!is.na(suppressWarnings(as.integer(param[i.par])))){
      names[i.par] = names(wlm.mdl$coefficients)[as.integer(param[i.par])]
    }
  }
  
  fitted.coef = coef(wlm.mdl)[names]
  # infl.vals computes w_i * psi(D_i)
  infl.vals = matrix(n * influence(wlm.mdl)$coefficients[,names], ncol=length(param))
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
      fit.Zr.new = (predict(Zr.1, data.frame(newdata)) + predict(Zr.2, data.frame(newdata)))/2
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
      fit.Zr.new = (predict(Zr.1, data.frame(newdata))$predictions + 
                      predict(Zr.2, data.frame(newdata))$predictions)/2
    }
    
    ##======================================##
    ##=== compute transductive estimator ===##  
    trans.coefs[i.par] = fitted.coef[i.par] - mean(ws * fit.Zr, na.rm=TRUE) + mean(fit.Zr.new, na.rm=TRUE)
    
    ##=======================================##
    ##== compute super-population variance ==##   
    sup.sds[i.par] = sqrt(sd(ws * (uw.infl.vals[,i.par] - fit.Zr), na.rm=TRUE)^2 + n/m * sd(fit.Zr.new, na.rm=TRUE)^2)
    
    ##===============================================##
    ##== compute conditional-transductive variance ==##    
    if (is.null(cond.newdata)){
      trans.sds[i.par] = min(sup.sds[i.par], sd(ws * (uw.infl.vals[,i.par] - fit.Zr), na.rm=TRUE))
    }else if(ncol(cond.newdata) == ncol(Z)){
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
        Xr.1 = loess(fit.Zr.new[fold2[fold2>n]]~., data=data.frame(cond.newdata[fold2[fold2>n],]),
                     span=loess.span, degree=loess.deg) 
        Xr.2 = loess(fit.Zr.new[fold1[fold1>n]]~., data=data.frame(cond.newdata[fold1[fold1>n],]),
                     span=loess.span, degree=loess.deg) 
        fit.Xr[fold1[fold1>n]] = predict(Xr.1, data.frame(cond.newdata[fold1[fold1>n],]))
        fit.Xr[fold2[fold2>n]] = predict(Xr.2, data.frame(cond.newdata[fold2[fold2>n],])) 
      } 
      # random forest regression
      if (alg == 'grf'){
        # cross-fit nonparametric regression models
        Zr.1 = regression_forest(data.frame(cond.newdata[fold2[fold2>n],]), 
                                 fit.Zr.new[fold2[fold2>n]], num.threads=1)
        Zr.2 = regression_forest(data=data.frame(cond.newdata[fold1[fold1>n],]), 
                                 fit.Zr.new[fold1[fold1>n]], num.threads=1)
        fit.Xr[fold1[fold1>n]] = predict(Zr.1, data.frame(cond.newdata[fold1[fold1>n],]))$predictions
        fit.Xr[fold2[fold2>n]] = predict(Zr.2, data.frame(cond.newdata[fold2[fold2>n],]))$predictions 
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
  sup.p_vals = 2 * pnorm(abs(trans.coefs) / (sup.sds/sqrt(n)), lower.tail = FALSE)
  
  # summary of transductive model
  ret_table = cbind(trans.coefs, trans.sds, trans.p_vals, sup.sds, sup.p_vals)
  colnames(ret_table) = c("Trans. Estimate", #"W. OLS Estimate", 
                          "Trans. Std. Error", "Trans. Pr(>|z|)",
                          "Sup. Std. Error", "Sup. Pr(>|z|)")
  rownames(ret_table) = names
  
  # print the summary if verbose==TRUE
  if (verbose){
    cat("\n")
    if (class(object)[1]=='lm'){
      cat("Summary of transductive inference in linear models")
    }else if (class(object)[1]=='glm'){
      cat("Summary of transductive inference in generalized linear models")
    } 
    cat("\n\n")
    print(ret_table)
    cat("\n")
  }
  
  invisible(list("trans.coef" = trans.coefs, #"OLS.coef" = fitted.coef, 
                 "trans.std.err" = trans.sds, "trans.pval" = trans.p_vals,
                 "trans.ci.low" = trans.ci.lo, "trans.ci.upp" = trans.ci.hi,
                 "sup.std.err" = sup.sds, "sup.pval" = sup.p_vals, "summary" = ret_table))
  
}


# for preparation, not exported
lm.internal <- function(formula, data, weights)
{
  envir <- list2env(list(weights=weights), parent=environment(formula))
  environment(formula) <- envir
  lm(formula, data=data, weights=weights)
}

glm.internal <- function(formula, data, family, weights)
{
  envir <- list2env(list(weights=weights), parent=environment(formula))
  environment(formula) <- envir
  glm(formula, data=data, family=family, weights=weights)
}