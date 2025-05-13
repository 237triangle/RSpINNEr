#' Prox F via SVD
#'
#' @param y
#' @param SVDAx
#' @param Dk
#' @param Wk
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
ProxFsvd <- function(y, SVDAx, Dk, Wk, delta) {

  # Proximity operator for the function F
  # F(B) := sum_i (y_i - <A_i, B>)^2
  # where A is assumed to be a 3-way tensor and y is the vector of observations.

  p <- nrow(Wk)

  # Extract components from SVDAx
  AU <- SVDAx$U
  ASdiag <- SVDAx$Sdiag
  AVt <- SVDAx$Vt
  idxs <- SVDAx$idxs

  # Compute Ddelta and related vectors
  Ddelta <- Dk + Wk / delta
  DdeltaVec <- as.vector(Ddelta)
  DdeltaVecU <- DdeltaVec[idxs]
  Mdelta <- t(AU) %*% DdeltaVecU * ASdiag
  ySubs <- t(AVt) %*% Mdelta
  ydelta <- y - ySubs

  # Ridge solution from SVD
  AVty <- AVt %*% ydelta
  MAVty <- (AVty * ASdiag) / (ASdiag^2 + 2 * delta)  # element-wise
  BridgeVec <- AU %*% MAVty

  # Reconstruct the symmetric matrix
  Bridge <- rep(0, p^2)
  Bridge[idxs] <- BridgeVec
  Bridge <- matrix(Bridge, nrow = p, ncol = p)  # p-by-p matrix
  Bridge <- Bridge + t(Bridge)  # Symmetric part

  # B update
  Bnew <- Bridge + Ddelta

  return(Bnew)
}

#' Prox G
#'
#' @param D
#' @param W2
#' @param delta
#' @param lambda_N
#'
#' @return
#' @export
#'
#' @examples
ProxG <- function(D, W2, delta, lambda_N) {

  # Proximity operator for the function G
  # G(C) := lambda_N * ||C||_*  (nuclear norm of C)

  # Compute Ddelta and symmetrize it
  Ddelta <- D + W2 / delta
  Ddelta <- (Ddelta + t(Ddelta)) / 2  # Ensuring symmetry

  # Perform eigendecomposition
  eig_res <- eigen(Ddelta)
  U <- eig_res$vectors
  S <- eig_res$values

  # Soft thresholding of singular values
  S_thresh <- sign(S) * pmax(abs(S) - lambda_N / delta, 0)

  # Reconstruct the matrix using the thresholded singular values
  Cnew <- U %*% diag(S_thresh) %*% t(U)

  return(Cnew)
}

#' Prox H
#'
#' @param B
#' @param C
#' @param delta1
#' @param delta2
#' @param W1
#' @param W2
#' @param lambda_L
#' @param WGTs
#'
#' @return
#' @export
#'
#' @examples
ProxH <- function(B, C, delta1, delta2, W1, W2, lambda_L, WGTs) {
  deltas     <- delta1 + delta2
  Bdelta1    <- B - W1 / delta1
  Bdelta2    <- C - W2 / delta2
  Bdelta     <- (delta1 * Bdelta1 + delta2 * Bdelta2) / deltas
  Dnew       <- sign(Bdelta) * pmax(abs(Bdelta) - WGTs * lambda_L / deltas, 0)
  return(Dnew)
}

