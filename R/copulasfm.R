
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Main Function
copSFM=function(Y,X,family,RHO,LB,UB,verbose=FALSE){

  like<-function(theta,Y,X,family){


    midX=X
    midY=Y

    XX=as.matrix(cbind(1,midX))
    K=ncol(XX)
    sigmav=abs(theta[K+1])
    sigmau=abs(theta[K+2])
    rho=theta[K+3]
    n=length(midY)
    m=n


    w=c(midY-t(theta[1:K])%*%t(XX))
    u=replicate(n,rtruncnorm(n, a=0.001, b=Inf, mean = 0, sd = sigmau))
    W=t(replicate(n,w))
    gv=dnorm(u+W,mean=0,sd=sigmav)+0.000001
    gv=matrix(t(gv),nrow=n,ncol=n)
    Gv=pnorm(u+W,mean=0,sd=sigmav)

    Fu=ptruncnorm(u, a=0.0001, b=Inf, mean = 0, sd = sigmau)

    Fu=c(abs(Fu))
    Gvv=c(abs(Gv))
    mm=length(Fu)
    for ( i in 1:mm){
      if (is.infinite(Fu[i]))  # control for optimization
        Fu=0.0000000000000001
      if (is.infinite(Gvv[i]))  # control for optimization
        Gvv=0.000000000000001
      if (is.nan(Fu[i]))  # control for optimization
        Fu=0.00000000000001
      if (is.nan(Gvv[i]))  # control for optimization
        Gvv=0.00000000000001
    }

    if (family==2){
      rho2=theta[length(theta)]
      gaucopula=BiCopPDF(Fu, Gvv, family=family, par=rho, par2=rho2)+0.00000001
    }else{
      gaucopula=BiCopPDF(Fu, Gvv, family=family, par=rho, par2=0)+0.00000001}

    gaucopula=matrix(gaucopula,nrow=m,ncol=n)
    hw=sum(log(1/m*diag(gv%*%gaucopula)))
    if (is.infinite(hw))  # control for optimization
      hw<--n*100
    if (is.nan(hw))  # control for optimization
      hw<--n*100

    if (isTRUE(verbose)) {
      message(sprintf("Sum of log-likelihood for Copula-SFA: %.4f", hw))
    }


    return(hw) # log likelihood

  }

  ### End Function #############3
  #=================================================


  # start here
  # select familty  copula upper and lower bouubd ( look at CDVine package)
  family=family  # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....
  LB=LB      #lower bound
  UB=UB   # upper bound
  RHO=RHO   # any value in bound



  #Gaussian (-.99, .99)
  #Student t (-.99, .99)
  #Clayton (0.1, Inf)
  cc=coef(lm(Y~X))
  K=length(cc)

  if (family==2){
    lower =c(rep(-Inf,K),0.01,0.01,LB,4.1)
    upper =c(rep(Inf,K+2),UB,50)
    cc=rep(0.1,K)
    start0=c(cc,sigmav=1,sigmau=1,rho=RHO,df=4)
  }else{

    lower =c(rep(-Inf,K),0.01,0.01,LB)
    upper =c(rep(Inf,K+2),UB)
    start0=c(cc,sigmav=1,sigmau=1,rho=RHO)
  }


  model <- optim(start0,like,Y=Y,X=X,family=family,
                 control = list(maxit=10,fnscale=-1),method="L-BFGS-B",
                 lower =lower,upper =upper, hessian=TRUE )
  # table of results
  coef<- model$par
  k=length(coef)
  model$se <- sqrt(-diag(solve(model$hessian)))

  for(i in 1:k){
    if (is.nan(model$se[i]))  # control for optimization
      model$se[i] <- sqrt(-diag(solve(-model$hessian)))[i]
  }

  n=length(Y)
  S.E.= model$se
  (paramsWithTs = cbind (model$par , coef/S.E. ) )
  stat=coef/S.E.
  pvalue <- 2*(1 - pnorm(abs(stat)))
  result <- cbind(coef,S.E.,stat,pvalue)
  result
  BIC= -2*model$value+ (log(n)*length(coef))
  AIC = -2*model$value + 2*length(coef)


  output=list(
    result=result,
    AIC=AIC,
    BIC=BIC,
    Loglikelihood=model$value
  )
  output

}
