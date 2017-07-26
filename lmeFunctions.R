#bootLRT does a boostrap likelihood ratio test
#@param: model1 is the 'true' model (lmer object) under which data is simulated in the bootstrap procedure 
#@param: model2 is the competing model (lmer object)
#@param: nboot (integer) is the number of replicate bootstrap sample
#@output: the model outputs a list of 2 elements lRTboot (table of all the LRT statistics) and  P (the corresponding p value)
bootLRT=function(model1,model2,nboot){
  output=list()
  output$obsLRT=2*(logLik(model2)-logLik(model1))
  lRTboot=rep(NA,nboot)
  for(i in 1:nboot){
    yboot=(simulate(model1))
    lRTboot[i]=2*(logLik(refit(model2,newresp=yboot))-logLik(refit(model1,newresp=yboot)))
  }
  output$LRTboot=lRTboot
  output$P=mean(lRTboot>output$obsLRT,na.rm=T)
  return(output)
}

#BCconf returns CI using bias-corrected bootstrap method for a given parameter (see Efron and Tibishrani 1986)
#@param: thebootcol is a vector of bootstrap estimates for the focal parameter
#@param: thestat is the point estimate for the focal parameter
#@param: low and high are the CI limits, defaults is 95%
#@output: vector of CI limits
BCconf=function(thebootcol,thestat,low=0.025,high=0.975){
  nboot=length(thebootcol)#number of bootstrap replicates
  z0=qnorm(sum( thebootcol < thestat)/nboot)
  temp=sort(thebootcol)
  zalpha=qnorm(low)
  tt=pnorm(2*z0+zalpha)
  ooo=trunc(tt*nboot)
  confpoint1=temp[ooo]
  zalpha=qnorm(high)
  tt=pnorm(2*z0+zalpha)
  ooo=trunc(tt*nboot)
  confpoint2=temp[ooo+1]
  return(c(confpoint1,  confpoint2))
}

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}