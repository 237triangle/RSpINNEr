#' SpINNEr
#'
#' @param y
#' @param X
#' @param AA
#' @param lambdaN
#' @param lambdaL
#' @param W
#'
#' @return
#' @export
#'
#' @examples
spinner <- function(y, X, AA, lambdaN, lambdaL, W = NULL) {

  # This function solves the optimization problem using ADMM.
  #--------------------------------------------------------------------------
  # argmin_{B, beta} { 0.5 * sum_i ( y_i - X %*% beta - <A_i, B> )^2 +
  #                    lambda_N * || B ||_* + lambda_L * || vec(B o W) ||_1 }
  #--------------------------------------------------------------------------

  # Check symmetric matrices AA
  if (!all.equal(AA, aperm(AA, c(2, 1, 3)))) {
    stop("Matrices A_i`s must be symmetric")
  }

  # Check zeros on diagonals of AA
  # OLD
  # if (sum(diag(abs(apply(AA, 3, sum)))) > 0) {
  #   stop("Matrices A_i`s must have zeros on diagonals")
  # }
  # NEW
  if(sum(diag(abs(apply(AA,c(1,2), sum)))) > 0){
    stop("Maatrices A_i's must have zeros on diagonals")
  }


  # Check if X is provided
  if (is.null(X) || all(X == 0)) {
    XtXXt <- matrix(0, 1, length(y))
    X <- matrix(0, length(y), 1)
  } else {
    XtXXt <- solve(t(X) %*% X, t(X))
  }

  # Dimension checks
  if (dim(AA)[3] != length(y)) {
    stop("The third dimension of 3-way tensor containing A_i`s should be the same as the length of y")
  }
  if (nrow(X) != length(y)) {
    stop("Number of rows in X and the length of y should coincide")
  }

  p <- dim(AA)[1]
  n <- length(y)

  # Default W if not provided
  if (is.null(W)) {
    if (lambdaL > 0) {
      W <- matrix(1, p, p) - diag(1, p)
    } else {
      W <- matrix(1, p, p)
    }
  }

  # Get rid of X from the optimization problem
  H <- diag(n) - X %*% XtXXt
  AAmatrix <- matrix(NA, p^2, n)
  for (i in 1:n) {
    AAmatrix[, i] <- as.vector(AA[, , i])
  }
  # NEW AARON
  AAmatrix <- t(AAmatrix)
  AAtilde <- H %*% AAmatrix
  AAtilde <- t(AAtilde)
  AAtilde <- array(AAtilde, dim = c(p, p, n)) # X regressed out
  ytilde <- H %*% y # X regressed out

  # Solver options
  solOptions <- list(
    deltaInitial1 = 100,   # Initial step length for nuclear norm update
    deltaInitial2 = 100,   # Initial step length for LASSO update
    scaleStep = 1,         # Initial scale for delta updates
    ratioStep = 1,         # Initial ratio for delta updates
    mu = 10,               # Max acceptable ratio between convergence rates
    deltaInc = 2,          # Delta increase factor
    deltaDecr = 2,         # Delta decrease factor
    ratioInc = 2,          # Ratio increase factor
    ratioDecr = 2,         # Ratio decrease factor
    maxIters = 50000,      # Max iterations
    epsPri = 1e-6,         # Primary residual convergence tolerance
    epsDual = 1e-6         # Dual residual convergence tolerance
  )

  # SVD
  Avecs <- matrix(NA, p^2, n)
  for (i in 1:n) {
    Avecs[, i] <- as.vector(AAtilde[, , i])
  }
  idxs <- matrix(upper.tri(matrix(1, p, p)), nrow = p^2, ncol = 1)
  AvecsUp <- 2 * Avecs[idxs, , drop = FALSE] # (p^2-p)/2-by-n

  # Economy-size SVD
  svd_result <- svd(AvecsUp) # U is (p^2-p)/2-by-n, S is n-by-n, V is n-by-n
  U <- svd_result$u
  Sdiag <- svd_result$d
  Vt <- t(svd_result$v)

  # SVD objects
  SVDAx <- list(U = U, Sdiag = Sdiag, Vt = Vt, idxs = idxs)

  # Solver selection
  solverType <- as.integer(lambdaN > 0) + 2 * as.integer(lambdaL > 0) + 1

  # Solver
  switch(solverType,
         `1` = { out <- minNormEstim(ytilde, SVDAx) },
         `2` = { out <- spinnerNuclear(ytilde, SVDAx, lambdaN, solOptions) },
         `3` = { out <- spinnerLasso(ytilde, SVDAx, lambdaL, W, solOptions) },
         `4` = { out <- spinnerBoth(ytilde, SVDAx, lambdaN, lambdaL, W, solOptions) }
  )

  # Estimates
  estim <- out$B
  beta <- XtXXt %*% (y - AAmatrix %*% as.vector(estim))

  # Optimal value
  DlastVec <- as.vector(estim)
  DlastVecU <- DlastVec[idxs]
  MDlast <- crossprod(U, DlastVecU) * Sdiag
  optVal <- 0.5 * sum((ytilde - Vt %*% MDlast)^2) + lambdaN * sum(svd(estim)$d) + lambdaL * sum(abs(estim))

  # Output list
  out$optVal <- optVal
  out$beta <- beta

  return(out)
}



#' Spinner Nuclear Norm only, internal
#'
#' @param y
#' @param SVDAx
#' @param lambda_N
#' @param solOptions
#'
#' @return
#' @export
#'
#' @examples
spinnerNuclear <- function(y, SVDAx, lambda_N, solOptions) {

  # This function solves the optimization problem:
  # argmin_B {  0.5 * sum_i (y_i - <A_i, B>)^2 + lambda_N || B ||_*  }

  # Extract parameters
  p0 <- nrow(SVDAx$U)
  p <- (1 + sqrt(1 + 8 * p0)) / 2

  # Solver options
  deltaInitial <- solOptions$deltaInitial1
  mu <- solOptions$mu
  deltaInc <- solOptions$deltaInc
  deltaDecr <- solOptions$deltaDecr
  maxIters <- solOptions$maxIters
  epsPri <- solOptions$epsPri
  epsDual <- solOptions$epsDual

  # SVD components
  U <- SVDAx$U
  Sdiag <- SVDAx$Sdiag
  Vt <- SVDAx$Vt
  idxs <- SVDAx$idxs

  # Initial primal and dual matrix
  Ck <- matrix(0, p, p)
  Wk <- matrix(0, p, p)

  # ADMM loop
  delta <- deltaInitial
  counterr <- 0
  stop <- FALSE

  while (!stop) {
    # Solve using Proximal Operator functions
    Bnew <- ProxFsvd(y, SVDAx, Ck, Wk, delta)
    Cnew <- ProxG(Bnew, -Wk, delta, lambda_N)

    Wk <- Wk + delta * (Cnew - Bnew)
    rk <- Cnew - Bnew
    sk <- Cnew - Ck

    rknorm <- norm(rk, type = "F")
    Bnorm <- norm(Bnew, type = "F")
    rknormR <- rknorm / Bnorm
    sknorm <- norm(sk, type = "F")
    sknormR <- sknorm / norm(Ck, type = "F")

    Ck <- Cnew
    counterr <- counterr + 1

    # Stopping criteria
    if (counterr > 10) {
      if (rknorm > mu * sknorm) {
        delta <- deltaInc * delta
      } else {
        if (sknorm > mu * rknorm) {
          delta <- delta / deltaDecr
        }
      }
    }

    if (rknormR < epsPri && sknormR < epsDual) {
      stop <- TRUE
    }

    if (counterr > maxIters) {
      stop <- TRUE
    }

    if (Bnorm < 1e-16) {
      stop <- TRUE
      Bnew <- matrix(0, p, p)
      Cnew <- matrix(0, p, p)
    }
  }

  # Optimal value
  ClastVec <- as.vector(Cnew)
  ClastVecU <- ClastVec[idxs]
  MClast <- t(U) %*% ClastVecU * Sdiag
  optVal <- 0.5 * norm(y - Vt %*% MClast, type = "F")^2 + lambda_N * sum(svd(Cnew)$d)

  # Output list
  out <- list(
    optVal = optVal,
    count = counterr,
    delta = delta,
    Blast = Bnew,
    Clast = Cnew,
    B = Cnew
  )

  return(out)
}



#' Spinner Lasso only, internal
#'
#' @param y
#' @param SVDAx
#' @param lambdaL
#' @param WGTs
#' @param solOptions
#'
#' @return
#' @export
#'
#' @examples
spinnerLasso <- function(y, SVDAx, lambdaL, WGTs, solOptions) {
  # Extract problem dimensions
  p0 <- nrow(SVDAx$U)
  p <- (1 + sqrt(1 + 8 * p0)) / 2

  ## Solver options
  delta <- solOptions$deltaInitial2
  mu <- solOptions$mu
  deltaInc <- solOptions$deltaInc
  deltaDecr <- solOptions$deltaDecr
  maxIters <- solOptions$maxIters
  epsPri <- solOptions$epsPri
  epsDual <- solOptions$epsDual

  ## SVD components
  U <- SVDAx$U
  Sdiag <- SVDAx$Sdiag
  Vt <- SVDAx$Vt
  idxs <- SVDAx$idxs

  ## Initialization
  Dk <- matrix(0, p, p)
  Wk <- matrix(0, p, p)

  ## ADMM loop
  counter <- 0
  stop <- FALSE

  while (!stop) {
    Bnew <- ProxFsvd(y, SVDAx, Dk, Wk, delta)
    Dnew <- ProxH_lasso(Bnew, delta, Wk, lambdaL, WGTs)
    Wk <- Wk + delta * (Dnew - Bnew)

    rk <- Dnew - Bnew
    sk <- Dnew - Dk
    rknorm <- norm(rk, type = "F")
    Bnorm <- norm(Bnew, type = "F")
    rknormR <- rknorm / Bnorm
    sknorm <- norm(sk, type = "F")
    sknormR <- sknorm / norm(Dk, type = "F")

    Dk <- Dnew
    counter <- counter + 1

    ## Stopping criteria
    if (counter %% 10 == 0) {
      if (rknorm > mu * sknorm) {
        delta <- delta * deltaInc
      } else if (sknorm > mu * rknorm) {
        delta <- delta / deltaDecr
      }
    }

    if (rknormR < epsPri && sknormR < epsDual) {
      stop <- TRUE
    }
    if (counter > maxIters) {
      stop <- TRUE
    }
    if (Bnorm < 1e-16) {
      stop <- TRUE
      Bnew <- matrix(0, p, p)
      Dnew <- matrix(0, p, p)
    }
  }

  ## Compute optimal value
  DlastVec <- as.vector(Dnew)
  DlastVecU <- DlastVec[idxs]
  MDlast <- crossprod(U, DlastVecU) * Sdiag
  residual <- y - crossprod(Vt, MDlast)
  optVal <- 0.5 * sum(residual^2) + lambdaL * sum(abs(Dnew))

  ## Output
  out <- list(
    optVal = optVal,
    count = counter,
    delta = delta,
    Blast = Bnew,
    Dlast = Dnew,
    B = Dnew
  )

  return(out)
}




#' Spinner without either penalty, internal
#'
#' @param y
#' @param SVDAx
#'
#' @return
#' @export
#'
#' @examples
minNormEstim <- function(y, SVDAx) {
  ## Problem size
  p0 <- nrow(SVDAx$U)
  p <- (1 + sqrt(1 + 8 * p0)) / 2

  ## Extract SVD components
  U <- SVDAx$U
  Sdiag <- SVDAx$Sdiag
  Vt <- SVDAx$Vt
  idxs <- SVDAx$idxs

  ## Estimate
  Vty <- (Vt %*% y) / Sdiag
  Bvec <- U %*% Vty

  Bestim <- numeric(p^2)
  Bestim[idxs] <- Bvec
  Bestim <- matrix(Bestim, nrow = p, ncol = p)
  Bestim <- Bestim + t(Bestim)

  ## Output
  out <- list(B = Bestim)
  return(out)
}


#' SpINNEr both, internal
#'
#' @param y
#' @param SVDAx
#' @param lambdaN
#' @param lambdaL
#' @param WGTs
#' @param solOptions
#'
#' @return
#' @export
#'
#' @examples
spinnerBoth <- function(y, SVDAx, lambdaN, lambdaL, WGTs, solOptions) {

  # This function solves the optimization problem using ADMM

  # Extract dimensions
  p0 <- nrow(SVDAx$U)
  p <- (1 + sqrt(1 + 8 * p0)) / 2

  # Solver options
  deltaInitial1 <- solOptions$deltaInitial1
  deltaInitial2 <- solOptions$deltaInitial2
  scaleStep <- solOptions$scaleStep
  ratioStep <- solOptions$ratioStep
  mu <- solOptions$mu
  deltaInc <- solOptions$deltaInc
  deltaDecr <- solOptions$deltaDecr
  ratioInc <- solOptions$ratioInc
  ratioDecr <- solOptions$ratioDecr
  maxIters <- solOptions$maxIters
  epsPri <- solOptions$epsPri
  epsDual <- solOptions$epsDual

  # SVD components
  U <- SVDAx$U
  Sdiag <- SVDAx$Sdiag
  Vt <- SVDAx$Vt
  idxs <- SVDAx$idxs

  # Initial primal and dual matrices
  Dk <- matrix(0, p, p)
  W1k <- matrix(0, p, p)
  W2k <- matrix(0, p, p)

  # ADMM loop
  delta1 <- deltaInitial1
  delta2 <- deltaInitial2
  counterr <- 0
  stop <- 0
  CsB <- numeric()
  DsB <- numeric()
  DsDp <- numeric()
  Dlts1 <- numeric()
  Dlts2 <- numeric()

  while (stop == 0) {
    # Proximity operators
    Bnew <- ProxFsvd(y, SVDAx, Dk, W1k, delta1)
    Cnew <- ProxG(Dk, W2k, delta2, lambdaN)
    Dnew <- ProxH(Bnew, Cnew, delta1, delta2, W1k, W2k, lambdaL, WGTs)

    # Update dual variables
    W1k <- W1k + delta1 * (Dnew - Bnew)
    W2k <- W2k + delta2 * (Dnew - Cnew)

    # Residuals
    rk1 <- Cnew - Bnew
    rk2 <- Dnew - Bnew
    sk <- Dnew - Dk
    rknorm1 <- norm(rk1, "F")
    Bnorm <- norm(Bnew, "F")
    rknormR1 <- rknorm1 / Bnorm
    rknorm2 <- norm(rk2, "F")
    rknormR2 <- rknorm2 / Bnorm
    sknorm <- norm(sk, "F")
    sknormR <- sknorm / norm(Dk, "F")
    counterr <- counterr + 1
    CsB[counterr] <- rknormR1
    DsB[counterr] <- rknormR2
    DsDp[counterr] <- sknormR
    Dlts1[counterr] <- delta1
    Dlts2[counterr] <- delta2
    Dk <- Dnew

    # Ratio update
    if (counterr %% 20 == 10) {
      if (rknorm1 > mu * rknorm2) {
        ratioStep <- ratioStep * ratioInc
      } else if (rknorm2 > mu * rknorm1) {
        ratioStep <- ratioStep / ratioDecr
      }
    }

    # Scale update
    if (counterr %% 20 == 0) {
      if (mean(c(rknorm1, rknorm2)) > mu * sknorm) {
        scaleStep <- scaleStep * deltaInc
      } else if (sknorm > mu * mean(c(rknorm1, rknorm2))) {
        scaleStep <- scaleStep / deltaDecr
      }
    }

    # Update deltas
    delta1 <- scaleStep * deltaInitial1
    delta2 <- scaleStep * ratioStep * deltaInitial2

    # Stopping criteria
    if (rknormR1 < epsPri && rknormR2 < epsPri && sknormR < epsDual) {
      stop <- 1
    }
    if (counterr > maxIters) {
      stop <- 1
    }
    if (Bnorm < 1e-16) {
      stop <- 1
      Bnew <- matrix(0, p, p)
      Cnew <- matrix(0, p, p)
      Dnew <- matrix(0, p, p)
    }
  }

  # Optimal value calculation
  DlastVec <- as.vector(Dnew)
  DlastVecU <- DlastVec[idxs]
  MDlast <- (t(U) %*% DlastVecU) * Sdiag
  optVal <- 0.5 * norm(y - t(Vt) %*% MDlast, "F")^2 + lambdaN * sum(svd(Dnew)$d) + lambdaL * sum(abs(Dnew))

  # Outputs
  out <- list()
  out$optVal <- optVal
  out$count <- counterr
  out$Dlts1 <- Dlts1
  out$Dlts2 <- Dlts2
  out$Blast <- Bnew
  out$Clast <- Cnew
  out$Dlast <- Dnew
  out$B <- Dnew

  return(out)
}
