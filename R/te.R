
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
## Compute technical efficiency
#
#EX: te1=TE1(coef,Y,X,family=family)


#' @export
TE1 <- function(theta, Y, X, family,
                nSim = 200,
                rho2 = NULL,
                seed = NULL,
                plot = FALSE) {

  if (!is.null(seed)) set.seed(seed)

  Y  <- as.numeric(Y)
  X  <- as.matrix(X)
  XX <- cbind(1, X)

  n  <- length(Y)
  K  <- ncol(XX)

  beta   <- theta[1:K]
  sigmav <- abs(theta[K + 1])
  sigmau <- abs(theta[K + 2])
  rho    <- theta[K + 3]

  if (is.null(rho2) && length(theta) >= (K + 4))
    rho2 <- theta[K + 4]

  # residual
  w <- as.vector(Y - XX %*% beta)

  # ----------------------------
  # Simulation draw: u ~ truncated normal
  # ----------------------------
  u_std <- matrix(
    truncnorm::rtruncnorm(
      n * nSim,
      a = 0,
      b = Inf,
      mean = 0,
      sd = 1
    ),
    nrow = n,
    ncol = nSim
  )

  u <- sigmau * u_std
  W <- matrix(w, n, nSim)

  # ----------------------------
  # Normal noise density
  # ----------------------------
  gv <- stats::dnorm(u + W, mean = 0, sd = sigmav)
  Gv <- stats::pnorm(u + W, mean = 0, sd = sigmav)

  Fu <- truncnorm::ptruncnorm(
    u,
    a = 0,
    b = Inf,
    mean = 0,
    sd = sigmau
  )

  gv[!is.finite(gv) | gv <= 0] <- 1e-12
  Gv[!is.finite(Gv)] <- 1e-12
  Fu[!is.finite(Fu)] <- 1e-12

  # ----------------------------
  # Copula density
  # ----------------------------
  cop <- .copula_pdf(
    u1 = as.vector(Fu),
    u2 = as.vector(Gv),
    family = family,
    rho = rho,
    rho2 = rho2
  )

  cop[!is.finite(cop) | cop <= 0] <- 1e-12
  cop_mat <- matrix(cop, n, nSim)

  wgt <- gv * cop_mat

  # ----------------------------
  # TE = E[exp(-u) | epsilon]
  # ----------------------------
  num <- rowMeans(exp(-u) * wgt)
  den <- rowMeans(wgt)
  den[den <= 1e-300] <- 1e-300

  te <- num / den

  # ----------------------------
  # Optional plot
  # ----------------------------
  if (plot) {
    graphics::plot(
      sort(te, decreasing = TRUE),
      type = "l",
      xlab = "Observation",
      ylab = "Technical Efficiency",
      main = "Technical Efficiency"
    )
  }

  return(te)
}
