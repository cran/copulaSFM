\name{sfa.simu}
\alias{sfa.simu}
\title{Simulate Data for Stochastic Frontier Analysis}
\usage{
sfa.simu(nob, alpha, sigV, sigU, family, rho, rho2,seed = NULL )
}

\arguments{
  \item{nob}{Number of observations to be simulated.}
  \item{alpha}{Vector of frontier coefficients including the intercept and slope parameters.}
  \item{sigV}{Standard deviation of noise term V.}
  \item{sigU}{Standard deviation of inefficiency term U.}
  \item{family}{Integer code specifying the copula family used to model the dependence between the noise term and the inefficiency term.}
  \item{rho}{Primary dependence parameter of the copula.}
  \item{rho2}{
    Optional second copula parameter (passed to \code{par2} in \code{VineCopula}).
    This must be provided for copula families that require a second parameter
    (e.g., degrees of freedom for the Student-t copula, or a shape parameter for BB copulas).
  }
    \item{seed}{
    Optional integer seed for reproducibility. If \code{NULL}, no seed is set.
  }
}

\value{
A list containing simulated output and inputs.
}

\description{
generates simulated data from a copula-based stochastic frontier model (SFM).
The function creates dependent error components using a specified copula family and then constructs output data according to the stochastic frontier specification.
}

\details{
The simulated data follow the stochastic frontier model

\deqn{
y_i = x_i' \alpha + V_i - U_i
}

where

\eqn{V_i \sim N(0,\sigma_V^2)} is a symmetric noise term and
\eqn{U_i \sim |N(0,\sigma_U^2)|} is a non-negative inefficiency term.

The dependence between \eqn{V_i} and \eqn{U_i} is modeled using a copula.
Let \eqn{(U_1,U_2)} be a pair of uniform random variables generated from a
copula \eqn{C_{\rho}(u_1,u_2)} with dependence parameter \eqn{\rho}.
The noise and inefficiency terms are obtained through the transformations

\deqn{
V_i = \Phi^{-1}(U_1)\sigma_V
}

\deqn{
U_i = F^{-1}_{TN}(U_2;0,\sigma_U^2)
}

where \eqn{\Phi^{-1}(\cdot)} is the inverse standard normal distribution
function and \eqn{F^{-1}_{TN}(\cdot)} is the inverse cumulative
distribution function of the truncated normal distribution with support
\eqn{[0,\infty)}.

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

\examples{
\donttest{
# Example: Simulate data from a copula-based stochastic frontier model
# Set a seed for reproducibility
set.seed(1)
sim <- sfa.simu(nob = 20, alpha = c(1, 0.5, -0.2), sigV = 1, sigU = 1, family = 1, rho = 0.2)
}
}
\keyword{models}
