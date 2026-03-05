
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#

# ------------------------------------------------------------
# Safe Copula Simulation (rho / rho2 interface)
# ------------------------------------------------------------
#' @noRd
.copula_sim <- function(n, family, rho, rho2 = NULL) {

  if (!is.null(rho2)) {
    chk <- VineCopula::BiCopCheck(family = family, par = rho, par2 = rho2)
    if (isTRUE(chk)) {
      return(
        VineCopula::BiCopSim(
          N = n,
          family = family,
          par = rho,
          par2 = rho2
        )
      )
    }
  }

  chk1 <- VineCopula::BiCopCheck(family = family, par = rho, par2 = 0)

  if (!isTRUE(chk1)) {
    stop("Invalid copula parameters.")
  }

  VineCopula::BiCopSim(
    N = n,
    family = family,
    par = rho,
    par2 = 0
  )
}

############## Simulation

sfa.simu <- function(nob,
                     alpha,
                     sigV,
                     sigU,
                     family,
                     rho,
                     rho2 = NULL,
                     seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  .sfa_internal(
    n = nob,
    alpha = alpha,
    sigV = sigV,
    sigU = sigU,
    family = family,
    rho = rho,
    rho2 = rho2
  )
}

.sfa_internal <- function(n, alpha, sigV, sigU, family, rho, rho2 = NULL) {

  U_sim <- .copula_sim(n = n, family = family, rho = rho, rho2 = rho2)

  V <- qnorm(U_sim[,1], mean = 0, sd = sigV)

  U_eff <- truncnorm::qtruncnorm(
    U_sim[,2],
    mean = 0,
    sd = sigU,
    a = 0,
    b = Inf
  )

  x1 <- rnorm(n,0,1)
  x2 <- rnorm(n,0,1)

  XX <- as.matrix(cbind(1,x1,x2))

  y <- c(t(alpha) %*% t(XX) + V - U_eff)

  XX <- as.matrix(cbind(x1,x2))

  out <- list(
    Y = y,
    X = XX
  )

  return(out)
}

