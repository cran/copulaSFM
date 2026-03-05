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
TE1(theta, Y, X, family,
                nSim = 200,
                rho2 = NULL,
                seed = NULL,
                plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{theta}{
    Numeric vector of estimated model parameters. The expected ordering is:
    regression coefficients (including intercept), \eqn{\sigma_v}, \eqn{\sigma_u},
    copula parameter \code{rho}, and optionally \code{rho2} if a two-parameter copula is used.
  }

  \item{Y}{
    Numeric vector of dependent variable observations.
  }

  \item{X}{
    Numeric matrix (or object coercible to a matrix) of independent variables
    (without intercept column).
  }

  \item{family}{
    Integer specifying the copula family (see \code{VineCopula::BiCopPDF}).
    Some families require a second parameter \code{rho2}.
  }

   \item{nSim}{
    Number of Monte Carlo draws used to approximate the conditional expectation.
    Larger values reduce simulation noise but increase computation time.
  }

  \item{rho2}{
    Optional second copula parameter. If \code{NULL} and \code{theta} contains an additional
    element beyond \code{rho}, the function will use that element as \code{rho2}.
  }

  \item{seed}{
    Optional integer seed for reproducibility in simulation-based computation.
    If \code{NULL}, no seed is set.
  }

  \item{plot}{
    Logical. If \code{TRUE}, produces a plot of sorted technical efficiency values.
  }


}
\details{

Technical efficiency is computed as a simulated conditional expectation:
\deqn{TE_i = E[\exp(-U_i) \mid w_i],}
where \eqn{w_i = Y_i - x_i^\beta} and the weights are constructed using the
asymmetric Laplace density for the noise term and the copula density capturing dependence
between the inefficiency component and the noise component.

The inefficiency term \eqn{U} is generated from a truncated asymmetric Laplace distribution
on \eqn{(0,\infty)}.

If \code{plot = TRUE}, the function produces a line plot of sorted technical efficiency values.

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
\describe{
A numeric vector of technical efficiency values (one per observation).
}


}
\references{
Wiboonpongse, A., Liu, J., Sriboonchitta, S., & Denoeux, T. (2015). Modeling dependence between error components of the stochastic frontier model using copula: application to intercrop coffee production in Northern Thailand. International Journal of Approximate Reasoning, 65, 34-44.

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
data=sfa.simu(nob=50, alpha=c(1,2,0.5),sigV=1,sigU=0.5,family=1,rho=0.5)

# Estimate the copula-based stochastic frontier model
model=copSFM(Y=data$Y,X=data$X,family=1,RHO=0.5,LB=-0.99,UB=0.99,  nSim = 50,  seed = 123)

## -------------------------------------------------
## Example: Computing Technical Efficiency
## -------------------------------------------------

# Compute technical efficiency estimates
te1=TE1(model$result[,1],Y=data$Y,X=data$X,family=1, nSim = 50,  seed = 123,  plot = TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{technical efficiency}
