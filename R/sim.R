
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
############## Simulation

sfa.simu<-function(nob,alpha,sigV,sigU,family,rho)
{
  nmax<-nob

  if (family==2){
    sim = VineCopula::BiCopSim(nob,family=family,rho,4)
  }else{
    sim = VineCopula::BiCopSim(nob,family=family,rho)
  }

  W=truncnorm::qtruncnorm(sim[,2], a=0.001, b=Inf, mean = 0, sd = sigU)

  V=(qnorm(sim[,1], mean=0,sd=sigV))

  # e <- rnorm(nob, 0, (sig))
  x1=rnorm(nmax,0,1)
  x2=rnorm(nmax,0,1)
  XX=as.matrix(cbind(1,x1,x2))
  y =c(t(alpha)%*%t(XX)+V-W)

  XX=as.matrix(cbind(x1,x2))
  out=list(Y=y,X=XX)
  return(out)
}

