\name{copSFM}
\alias{copSFM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copula based Stochastic frontier Model
}
\description{
In the standard stochastic frontier model, the two-sided error term V and the one-sided technical inefficiency error term W are assumed to be independent. In this paper, we relax this assumption by modeling the dependence between V and W using copulas.
}

\usage{
copSFM(Y, X, family,
                   RHO, LB, UB,
                   RHO2 = NULL, LB2 = NULL, UB2 = NULL,
                   nSim = 50,
                   seed = NULL,
                   maxit = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Y}{ vector of dependent variable}
\item{X}{ matrix of independent variable}
\item{family}{ Copula function eg. Gaussain=1, Student-t=2 (see, Vinecopula package)}
\item{RHO}{ Initial value of the first copula parameter (passed to \code{par} in \code{VineCopula}).
    Note that this parameter is not always a correlation; its interpretation depends on
    the chosen copula family.}
\item{LB}{  The lower bound of the copula parameter}
\item{UB}{ The upper bound of the copula parameter}
\item{RHO2}{
    Optional initial value of the second copula parameter (passed to \code{par2} in \code{VineCopula}).
    This is required for copula families that have two parameters (e.g. Student-t where \code{RHO2}
    represents degrees of freedom, and BB families where \code{RHO2} is a shape parameter).
    If \code{NULL}, a one-parameter copula is assumed.
  }
  \item{LB2}{
    Optional lower bound for the second copula parameter \code{RHO2}.
    Must be supplied when \code{RHO2} is not \code{NULL}.
  }
  \item{UB2}{
    Optional upper bound for the second copula parameter \code{RHO2}.
    Must be supplied when \code{RHO2} is not \code{NULL}.
  }
  \item{nSim}{
    Number of Monte Carlo draws used to approximate the likelihood integral for each observation.
    Larger values improve accuracy but increase computation time.
  }
  \item{seed}{
    Optional integer seed for reproducibility of the simulation draws used in the likelihood.
    If \code{NULL}, no seed is set.
  }
  \item{maxit}{
    Maximum number of iterations for the optimizer (\code{stats::optim} with \code{"L-BFGS-B"}).
}
}

\details{
The model follows the stochastic frontier decomposition
\deqn{Y_i = x_i\beta + v_i - u_i,}
where \eqn{v_i} is a two-sided noise term and \eqn{u_i \ge 0} is the inefficiency term.
Dependence between \eqn{u_i} and \eqn{v_i} is introduced through a bivariate copula density
\eqn{c(\cdot,\cdot)} from the \code{VineCopula} package. The log-likelihood is evaluated by
simulating draws of \eqn{u_i} from the truncated ALD and approximating the likelihood integral by
Monte Carlo averaging (\code{nSim} draws per observation).

When a two-parameter copula family is used, \code{RHO2} must be provided along with bounds
\code{LB2} and \code{UB2}.

The following copula families are supported, together with their
parameter bounds:

\preformatted{
1  = Gaussian copula (par: (LB = -0.99, UB = 0.99))
2  = Student t copula (par: (LB = -0.99, UB = 0.99); par2: (LB = 0, UB = Inf))
3  = Clayton copula (par: (LB = 0.1, UB = Inf))
4  = Gumbel copula (par: [LB = 0.99, UB = Inf))
5  = Frank copula (par: (LB = -Inf, UB =  0) U (LB = 0, UB = Inf))
6  = Joe copula (par: (LB = 0.99, UB = Inf))
7  = BB1 (Clayton-Gumbel) copula (par: (LB = 0, UB = Inf); par2: [LB = 0.99, UB = Inf))
8  = BB6 copula (par: (LB = 0.99, UB = Inf); par2: (LB = 0, UB = Inf))
9  = BB7 (Joe-Clayton) copula (par: [LB = 0.99, UB = Inf); par2: (LB = 0, UB = Inf))
10 = BB8 copula (par: (LB = 0, UB = 1); par2: (LB = 0, UB =  Inf))
13 = Survival Clayton (180 degrees rotation; par: (LB = 0, UB = Inf))
14 = Survival Gumbel (180 degrees rotation; par: [LB = 0.99, UB = Inf))
16 = Survival Joe (180 degrees rotation; par: (LB = 0.99, UB = Inf))
17 = Survival BB1 (par: (LB = 0, UB = Inf); par2: [LB = 0.99, UB = Inf))
18 = Survival BB6 (par: (LB = 0.99, UB = Inf); par2: (0, UB = Inf))
19 = Survival BB7 (par: [LB = 0.99, UB = Inf); par2: (0, UB = Inf))
20 = Survival BB8 (par: (LB = 0, 0.99); par2: (LB = 0, UB = Inf))
23 = Rotated Clayton (90 degrees; par: (LB = -Inf, UB = 0))
24 = Rotated Gumbel (90 degrees; par: (LB = -Inf, UB = -0.99])
26 = Rotated Joe (90 degrees; par: (LB = -Inf, UB = -0.99))
27 = Rotated BB1 (90 degrees; par: (LB = -Inf, UB = 0); par2: [LB = 0.99, UB = Inf))
28 = Rotated BB6 (90 degrees; par: (LB = -Inf, UB = -0.99]; par2: (LB = 0, UB = Inf))
29 = Rotated BB7 (90 degrees; par: (LB = -Inf, UB = -0.99]; par2: (LB = 0, UB = Inf))
30 = Rotated BB8 (90 degrees; par: (LB = -0.99, UB = 0); par2: (LB = 0, UB = Inf))
33 = Rotated Clayton (270 degrees; par: (LB = -Inf, UB = 0))
34 = Rotated Gumbel (270 degrees; par: (LB = -Inf, UB = -0.99])
36 = Rotated Joe (270 degrees; par: (LB = -Inf, UB = -0.99))
37 = Rotated BB1 (270 degrees; par: (LB = -Inf, UB = 0); par2: [LB = 0.99, UB = Inf))
38 = Rotated BB6 (270 degrees; par: (LB = -Inf, UB = -0.99]; par2: (LB = 0, UB = Inf))
39 = Rotated BB7 (270 degrees; par: (LB = -Inf, UB = -0.99]; par2: (LB = 0, UB = Inf))
40 = Rotated BB8 (270 degrees; par: (LB = -0.99, UB = 0); par2: (LB = 0, UB = Inf))
}
See the \pkg{VineCopula} package (Nagler et al., 2025) for additional details on copula
parameterization.



}
\value{
\item{result }{The result contain the estimated parameters, standard errors, t-stat, and p-value}
\item{AIC }{Akaiki Information Criteria}
\item{BIC }{Bayesian Information Criteria}
\item{Loglikelihood }{Maximum Log-likelihood function}
\item{convergence}{
  Convergence code returned by \code{stats::optim}. A value of 0 indicates successful convergence.
}
}
\references{
Wiboonpongse, A., Liu, J., Sriboonchitta, S., & Denoeux, T.(2015). Modeling dependence between error components of the stochastic frontier model using copula: application to intercrop coffee production in Northern Thailand. International Journal of Approximate Reasoning, 65, 34-44.

Maneejuk, P., Yamaka, W., & Sriboonchitta, S.(2017). Analysis of global competitiveness using copula-based stochastic frontier kink model. In Robustness in Econometrics (pp. 543-559). Springer, Cham.

Nagler, T., Schepsmeier, U., Stoeber, J., Brechmann, E. C., Graeler, B., Erhardt, T., ... & Killiches, M. (2025, July). Package ‘VineCopula’.
}
\author{
Woraphon Yamaka, Paravee Maneejuk, and Nuttaphong Kaewtathip
}

\examples{
\donttest{
## -------------------------------------------------
## Example: Simulation and Model Estimation
## -------------------------------------------------

# Generate simulated stochastic frontier data
data=sfa.simu(nob=50, alpha=c(1,2,0.5),sigV=1,sigU=0.5,family=1,rho=0.5)

## -------------------------------------------------
## Example 1: One-parameter copula (Gaussian)
## -------------------------------------------------

# Specify the copula family and corresponding parameter bounds.
# For example:
# family = 1  -> Gaussian copula   (parameter range: -0.99, 0.99)
# family = 2  -> Student-t copula  (parameter range: -0.99, 0.99)
# family = 3  -> Clayton copula    (parameter range:  0.1,  Inf)

model=copSFM(Y=data$Y,X=data$X,family=1,RHO=0.5,LB=-0.99,UB=0.99)

## -------------------------------------------------
## Example 2: Two-parameter copula (Student-t)
## -------------------------------------------------

# For two-parameter copulas (e.g., Student-t, BB families),
# the second parameter (RHO2) must be supplied along with
# its corresponding lower and upper bounds.
# For the Student-t copula, RHO2 typically represents
# the degrees of freedom parameter.

data=sfa.simu(nob=50, alpha=c(1,2,0.5),sigV=1,sigU=0.5,family=2,rho=0.3,rho2=4)

model_t <- copSFM(
  Y = data$Y,
  X = data$X,
  family = 2,
  RHO = 0.3,
  LB = -0.99,
  UB = 0.99,
  RHO2 = 45,
  LB2 = 1,
  UB2 = Inf,
  nSim = 50,
  seed = NULL
)
model_t

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Copula }% use one of  RShowDoc("KEYWORDS")
\keyword{ Stochastic frointer }% __ONLY ONE__ keyword per line
