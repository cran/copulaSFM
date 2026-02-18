
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Compute technical efficiency
#
#EX: te1=TE1(coef,Y,X,family=family)


TE1<-function(theta,Y,X,family){

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

  hw=1/n*diag(gv%*%gaucopula) # A
  tu=sapply(u,mean)
  tcltcopula=sapply(gaucopula,mean)
  tgv=sapply(t(gv),mean)
  hw2=exp(-tu)*tcltcopula*tgv
  hw2=matrix(hw2,nrow=n,ncol=n)
  hw2=colMeans(hw2) # B
  # technical efficiency TE=A/B
  te=hw2/hw
  n=length(te)
  plot(te,lty=1,col="white",xlab = 'observation', ylab = 'TE value', main = "Technical Efficiency")
  lines(te, lty=1,type="l",col="blue")
  return(te)
}
