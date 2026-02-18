\name{sfa.simu}
\alias{sfa.simu}
\title{Simulate Data for Stochastic Frontier Analysis}
\usage{
sfa.simu(nob, alpha, sigV, sigU, family, rho)
}
\arguments{
  \item{nob}{Number of observations.}
  \item{alpha}{Coefficient vector.}
  \item{sigV}{Standard deviation of noise term V.}
  \item{sigU}{Standard deviation of inefficiency term U.}
  \item{family}{Copula family code.}
  \item{rho}{Copula dependence parameter.}
}
\value{
A list containing simulated output and inputs.
}
\description{
Simulates data for the copula-based stochastic frontier model.
}
\examples{
\donttest{
set.seed(1)
sim <- sfa.simu(nob = 20, alpha = c(1, 0.5, -0.2), sigV = 1, sigU = 1, family = 1, rho = 0.2)
}
}
\keyword{models}
