\name{copSFM}
\alias{copSFM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copula based Stochastic frontier Model
}
\description{
In the standard stochastic frontier model, the two-sided error term V and the one-sided technical inefficiency error term W are assumed to be independent. In this paper, we relax this assumption by modeling the dependence between V and W using copulas. Nine copula families are considered and their parameters are estimated using maximum simulated likelihood.
}
\usage{
copSFM(Y,X,family,RHO,LB,UB,verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y}{ vector of dependent variable}
\item{X}{ matrix of independent variable}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
\item{RHO}{ The initail value of the copula parameter}
\item{LB}{  The lower bound of the copula parameter}
\item{UB}{ The upper bound of the copula parameter}
\item{verbose}{ Logical; if TRUE, prints progress messages during optimization.}
}
\details{
herefore, the above copula families and relevant rotated copula can potentially capture the appropriate dependence between two random variables. Other popular copula families, such as Gaussain, Student,t Clayton, Gumbel etc.

}
\value{
\item{result }{The result contain the estimated parameters, standard errors, t-stat, and p-value}
\item{AIC }{Akaiki Information Criteria}
\item{BIC }{Bayesian Information Criteria}
\item{Loglikelihood }{Maximum Log-likelihood function}
}
\references{
Wiboonpongse, A., Liu, J., Sriboonchitta, S., & Denoeux, T.(2015). Modeling dependence between error components of the stochastic frontier model using copula: application to intercrop coffee production in Northern Thailand. International Journal of Approximate Reasoning, 65, 34-44.

Maneejuk, P., Yamaka, W., & Sriboonchitta, S.(2017). Analysis of global competitiveness using copula-based stochastic frontier kink model. In Robustness in Econometrics (pp. 543-559). Springer, Cham.
}
\author{
Woraphon Yamaka and Paravee MAneejuk
}

\examples{

#example simulation data
data=sfa.simu(nob=50, alpha=c(1,2,0.5),sigV=1,sigU=0.5,family=1,rho=0.5)

# Select familty  copula upper and lower bouubd ( look at CDVine package)
# family=1   # 1 is Gaussian, 2 is Student-t, 3 is Clayton and so on....

#Gaussian (-.99, .99)
#Student t (-.99, .99)
#Clayton (0.1, Inf)
model=copSFM(Y=data$Y,X=data$X,family=1,RHO=0.5,LB=-0.99,UB=0.99)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Stochastic frointer }% __ONLY ONE__ keyword per line
