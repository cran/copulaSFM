\name{TE1}
\alias{TE1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Technical efficiency measure.
}
\description{
Computing and plotting the technical efficiency.
}
\usage{
TE1(theta,Y,X,family)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{theta}{ The estimated parameters form the model}
 \item{Y}{ Vector of dependent variable}
\item{X}{ Matrix of independent variable}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
}
\details{
Computing and plotting the technical efficiency.
}
\value{
\describe{
  \item{Output}{Technical efficiency series.}
  \item{plot}{Plot of technical efficiency.}
}


}
\references{
Wiboonpongse, A., Liu, J., Sriboonchitta, S., & Denoeux, T. (2015). Modeling dependence between error components of the stochastic frontier model using copula: application to intercrop coffee production in Northern Thailand. International Journal of Approximate Reasoning, 65, 34-44.
}
\author{
Woraphon Yamaka
}

\examples{
## Required packages
#example simulation data
data=sfa.simu(nob=50, alpha=c(1,2,0.5),sigV=1,sigU=0.5,family=1,rho=0.5)

# Select familty  copula upper and lower bouubd ( look at CDVine package)
# family=1   # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....

#Gaussian (-.99, .99)
#Student t (-.99, .99)
#Clayton (0.1, Inf)

model=copSFM(Y=data$Y,X=data$X,family=1,RHO=0.5,LB=-0.99,UB=0.99)

#EX: Plot the technical efficiency
te1=TE1(model$result[,1],Y=data$Y,X=data$X,family=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
