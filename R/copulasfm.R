
# This is an  function named 'Copula based Stochastic frontier'
# Written by Dr.Worphon Yamaka,
# Center of excellence in Econometrics, Faculty of Economics,
# Chiang Mai University
#
#
#' @noRd
.copula_pdf <- function(u1, u2, family, rho, rho2 = NULL) {

  # families that REQUIRE par2
  two_par_fam <- c(2, 7, 8, 9, 10,
                   17,18,19,20,
                   27,28,29,30,
                   37,38,39,40)

  if (family %in% two_par_fam) {

    if (is.null(rho2)) {
      stop("This copula family requires two parameters (par2 missing).")
    }

    chk <- VineCopula::BiCopCheck(
      family = family,
      par    = rho,
      par2   = rho2
    )

    if (!isTRUE(chk)) {
      stop("Invalid copula parameters for this two-parameter family.")
    }

    return(
      VineCopula::BiCopPDF(
        u1, u2,
        family = family,
        par  = rho,
        par2 = rho2
      )
    )
  }

  # ----- one-parameter families -----

  chk1 <- VineCopula::BiCopCheck(
    family = family,
    par    = rho,
    par2   = 0
  )

  if (!isTRUE(chk1)) {
    stop("Invalid copula parameters.")
  }

  VineCopula::BiCopPDF(
    u1, u2,
    family = family,
    par  = rho,
    par2 = 0
  )
}




## Main Function
#' @importFrom stats coef
#' @export
copSFM <- function(Y, X, family,
                   RHO, LB, UB,
                   RHO2 = NULL, LB2 = NULL, UB2 = NULL,
                   nSim = 50,
                   seed = NULL,
                   maxit = 10000) {

  Y  <- as.numeric(Y)
  X  <- as.matrix(X)
  XX <- cbind(1, X)

  nObs <- length(Y)
  K    <- ncol(XX)

  if (!is.null(seed)) set.seed(seed)

  # ----------------------------
  # Simulation draw
  # ----------------------------
  u_std <- matrix(
    truncnorm::rtruncnorm(
      nObs * nSim,
      a = 0,
      b = Inf,
      mean = 0,
      sd = 1
    ),
    nrow = nObs,
    ncol = nSim
  )

  # ----------------------------
  # Log-likelihood
  # ----------------------------
  like <- function(theta) {

    beta   <- theta[1:K]
    sigmav <- abs(theta[K + 1])
    sigmau <- abs(theta[K + 2])
    rho    <- theta[K + 3]
    rho2 <- if (length(theta) > (K + 3)) theta[K + 4] else NULL

    w <- as.vector(Y - XX %*% beta)

    u <- sigmau * u_std
    W <- matrix(w, nObs, nSim)

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

    cop <- .copula_pdf(
      u1 = as.vector(Fu),
      u2 = as.vector(Gv),
      family = family,
      rho = rho,
      rho2 = rho2
    )

    cop[!is.finite(cop) | cop <= 0] <- 1e-12
    cop_mat <- matrix(cop, nObs, nSim)

    simLik <- rowMeans(gv * cop_mat)
    simLik[!is.finite(simLik) | simLik <= 0] <- 1e-12

    sum(log(simLik))
  }

  # ----------------------------
  # Starting values
  # ----------------------------
  cc <- stats::coef(stats::lm(Y ~ X))

  if (is.null(RHO2)) {

    lower  <- c(rep(-Inf, K), 0.01, 0.01, LB)
    upper  <- c(rep( Inf, K), Inf,  Inf,  UB)
    theta0 <- c(cc, 1, 1, RHO)

  } else {

    if (is.null(LB2) || is.null(UB2))
      stop("If RHO2 is used, LB2 and UB2 must be supplied.")

    lower  <- c(rep(-Inf, K), 0.01, 0.01, LB, LB2)
    upper  <- c(rep( Inf, K), Inf,  Inf,  UB, UB2)
    theta0 <- c(cc, 1, 1, RHO, RHO2)
  }

  # ----------------------------
  # Optimization
  # ----------------------------
  model <- stats::optim(
    theta0,
    like,
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(fnscale = -1, maxit = maxit),
    hessian = TRUE
  )

  coef_est <- model$par
  H        <- model$hessian

  vcov_mat <- tryCatch(
    solve(-H),
    error = function(e) MASS::ginv(-H)
  )

  se   <- sqrt(diag(vcov_mat))
  zval <- coef_est / se
  pval <- 2 * (1 - stats::pnorm(abs(zval)))

  # ----------------------------
  # Parameter names
  # ----------------------------
  beta_names <- colnames(X)

  if (is.null(beta_names))
    beta_names <- paste0("X", seq_len(ncol(X)))

  beta_names <- c("(Intercept)", beta_names)

  par_names <- c(
    beta_names,
    "sigma_v",
    "sigma_u",
    "rho"
  )

  if (!is.null(RHO2))
    par_names <- c(par_names, "rho2")

  result <- cbind(
    Estimate  = coef_est,
    Std.Error = se,
    z.value   = zval,
    p.value   = pval
  )

  rownames(result) <- par_names

  AIC <- -2 * model$value + 2 * length(coef_est)
  BIC <- -2 * model$value + log(nObs) * length(coef_est)

  output <- list(
    result        = result,
    AIC           = AIC,
    BIC           = BIC,
    Loglikelihood = model$value,
    convergence   = model$convergence
  )

  class(output) <- "copSFM"
  return(output)
}

### End Function #############3
#=================================================



