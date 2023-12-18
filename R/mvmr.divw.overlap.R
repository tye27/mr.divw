#'  Perform debiased adjusted inverse-variance weighted (adIVW) estimator for summary-data multivariable Mendelian randomization allowing for overlapping exposure and outcome datasets
#'
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param gen_cor Provide a (K+1)-by-(K+1) matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure and the last index position corresponds to the outcome. The correlations can either be estimated, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. Default input is NULL, meaning that an identity matrix is used as the correlation matrix.
#' @param phi_cand A vector of tuning parameters for adIVW estimator. Default is 0, meaning that dIVW estimator is performed. To use the recommended set for the tuning parameter, simply set phi_cand = NULL.
#'
#' @return A list with elements
#' \item{beta.hat}{Estimated direct effects of each exposure on the outcome}
#' \item{beta.se}{Estimated standard errors of beta.hat}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies the IV strength in the sample}
#' \item{phi_selected}{The selected tuning parameter for the adIVW estimator}
#' @import MVMR
#' @export
#'
#' @examples
#' library(MVMR)
#' data("rawdat_mvmr")
#' beta.exposure <- rawdat_mvmr[,c("LDL_beta","HDL_beta","Trg_beta")]
#' se.exposure <- rawdat_mvmr[,c("LDL_se","HDL_se","Trg_se")]
#' beta.outcome <- rawdat_mvmr$SBP_beta
#' se.outcome <- rawdat_mvmr$SBP_se
#' P <- matrix(0.3, nrow = 4, ncol = 4)
#' diag(P) <- 1
#' mvmr.divw(beta.exposure = beta.exposure,
#' se.exposure = se.exposure,
#' beta.outcome = beta.outcome,
#' se.outcome = se.outcome,
#' gen_cor = P,
#' phi_cand = NULL,
#' over.dispersion = FALSE)
#'
mvmr.divw.overlap <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = NULL, phi_cand=0) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {stop("this function is developed for multivariable MR")}
  K <- ncol(beta.exposure)
  if (is.null(gen_cor)) {
    P <- diag(K)
  } else {
    P <- as.matrix(gen_cor)
  }
  if (nrow(beta.exposure) != length(beta.outcome)) {stop("The number of SNPs in beta.exposure and beta.outcome is different")}
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  # the number of SNPs
  p <- nrow(beta.exposure)
  # diagonal W matrix
  W<- diag(se.outcome^(-2))
  # create a list of Sigma Xj matrices
  Vj <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]))
  # calculate square root inverse of P
  P_eigen <- eigen(P)
  P_root_inv <- P_eigen$vectors %*% diag(1/sqrt(P_eigen$values)) %*% t(P_eigen$vectors)
  Vj_root_inv <- lapply(1:p, function(j) {
    P_root_inv %*% diag(1/se.exposure[j,])
  })
  # calcualte IV strenght parameter
  IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
    beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
  iv_strength_parameter <- min(eigen(IV_strength_matrix/sqrt(p))$values)
  # get V matrix
  V <- Reduce("+",lapply(1:p, function(j) {Vj[[j]] * (se.outcome[j]^(-2))}))
  # get M matrix
  M <- t(beta.exposure)%*%W%*%beta.exposure
  # get M-V matrix
  MV <- M-V
  # perform eigen-decomposition on MV
  MV_eigvalues <- eigen(MV)$values
  MV_eigen <- eigen(MV)
  # multivariable (a)dIVW estimator
  if (is.null(phi_cand)) {
    phi_cand <- c(0, exp(seq(0, 15, by = 0.5) - min(iv_strength_parameter)/2))
  }
  phi_length <- length(phi_cand)
  MV.l.inv.long <- Reduce(rbind, lapply(1:phi_length, function(l) {
    MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_cand[l]/MV_eigvalues)) %*% t(MV_eigen$vectors)}
  ))
  # (a)dIVW allowing for overlapping datasets
  Adj_term <- Reduce("+",lapply(1:p, function(j) {Vj[[j]][1:K,(K+1)] * se.outcome[j]^{-2}}))
  beta.est <- MV.l.inv.long %*% t(beta.exposure) %*% W %*% (beta.outcome) - MV.l.inv.long %*% Adj_term
  prof.lik <- sapply(1:phi_length, function(l) {
    beta.hat <- beta.est[(1+(l-1)*K):(l*K)]
    bvb.test <- sapply(1:p, function(j) t(beta.hat) %*% Vj[[j]][1:K,1:K] %*% beta.hat)
    bvxy <- sapply(1:p, function(j) {
      sigmaxy <- Vj[[j]][1:K,K+1]
      (t(beta.hat) %*% sigmaxy)[1,1]
    })
    S <- diag(1/(se.outcome^2 + bvb.test - 2 * bvxy))
    1/p * t(beta.outcome - beta.exposure %*% beta.hat) %*% S %*%
      (beta.outcome - beta.exposure %*% beta.hat)})
  phi_selected <- phi_cand[which.min(prof.lik)]
  MV.l.inv <- MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_selected/MV_eigvalues)) %*% t(MV_eigen$vectors)
  mvmr.adIVW <- MV.l.inv %*% t(beta.exposure) %*% W %*% beta.outcome - MV.l.inv %*% Adj_term
  adIVW_Vt_overlap <- Reduce("+",lapply(1:p, function(j) {
    m <- beta.exposure[j,] %*% t(beta.exposure[j,]) * (se.outcome[j]^(-2))
    v <- Vj[[j]][1:K,1:K] * (se.outcome[j]^(-2))
    bvb <- as.numeric(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
    vbbv <- v %*% mvmr.adIVW %*% t(mvmr.adIVW) %*% v
    sigmaxy <- Vj[[j]][1:K,K+1]
    A1 <- sigmaxy %*% t(sigmaxy) * se.outcome[j]^(-4)
    A2 <- Vj[[j]][1:K,1:K] %*% mvmr.adIVW %*% t(sigmaxy) * se.outcome[j]^(-4)
    A3 <- sigmaxy %*% t(mvmr.adIVW) %*% Vj[[j]][1:K,1:K] * se.outcome[j]^(-4)
    A4 <- (t(mvmr.adIVW) %*% sigmaxy * se.outcome[j]^(-2))[1,1] *  m
    A5 <- (t(mvmr.adIVW) %*% sigmaxy * se.outcome[j]^(-2))[1,1] *  sigmaxy %*% t(sigmaxy) * se.outcome[j]^(-4)
    m*(1+bvb) + vbbv + A1 - A2 - A3 - 2 * A4 + 4 * A5
  }))
  mvmr.adIVW.se <- sqrt(diag(MV.l.inv%*%adIVW_Vt_overlap%*%MV.l.inv))
  return(list(beta.hat = mvmr.adIVW,
              beta.se = mvmr.adIVW.se,
              iv_strength_parameter = iv_strength_parameter,
              phi_selected = phi_selected,
              tau.square = tau2_adivw))
}
