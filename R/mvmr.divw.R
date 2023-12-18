#'  Perform debiased adjusted inverse-variance weighted (adIVW) estimator for two-sample summary-data multivariable Mendelian randomization
#'
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param gen_cor Provide a K-by-K matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure. The correlations can either be estimated, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. Default input is NULL, meaning that an identity matrix is used as the correlation matrix.
#' @param phi_cand A vector of tuning parameters for adIVW estimator. Default is 0, meaning that dIVW estimator is performed. To use the recommended set for the tuning parameter, simply set phi_cand = NULL.
#' @param over.dispersion Should the model consider balanced horizontal pleiotropy. Default is FALSE
#'
#' @return A list with elements
#' \item{beta.hat}{Estimated direct effects of each exposure on the outcome}
#' \item{beta.se}{Estimated standard errors of beta.hat}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies the IV strength in the sample}
#' \item{phi_selected}{The selected tuning parameter for the adIVW estimator}
#' \item{tau.square}{Overdispersion parameter if \code{over.dispersion=TRUE}}
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
#' P <- matrix(0.3, nrow = 3, ncol = 3)
#' diag(P) <- 1
#' mvmr.divw(beta.exposure = beta.exposure,
#' se.exposure = se.exposure,
#' beta.outcome = beta.outcome,
#' se.outcome = se.outcome,
#' gen_cor = P,
#' phi_cand = NULL,
#' over.dispersion = FALSE)
#'
mvmr.divw <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = NULL, phi_cand=0, over.dispersion = FALSE) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {stop("this function is developed for multivariable MR")}
  K <- ncol(beta.exposure)
  if (is.null(gen_cor)) {
    P <- diag(K)
  } else {
    P <- as.matrix(gen_cor)
  }
  if (ncol(P) != ncol(beta.exposure)) {stop("The shared correlation matrix has a different number of columns than the input beta.exposure")}
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
  # (a)dIVW for independent datasets
  beta.est <- MV.l.inv.long %*% t(beta.exposure) %*% W %*% (beta.outcome)
  prof.lik <- sapply(1:phi_length, function(l) {
        beta.hat <- beta.est[(1+(l-1)*K):(l*K)]
        bvb <- sapply(1:p, function(j) t(beta.hat) %*% Vj[[j]] %*% beta.hat)
        if (over.dispersion) {
          tau2 <- Reduce("+",lapply(1:p, function(j) {
                  v <- Vj[[j]] * (se.outcome[j]^(-2))
                  (beta.outcome[j] - beta.exposure[j,] %*% beta.hat)^2*se.outcome[j]^(-2) - 1 - as.numeric(t(beta.hat) %*% v %*% beta.hat)
        }))/sum(diag(W))
          tau2 <- as.numeric(tau2)
        } else {
          tau2 <- 0
        }
        S <- diag(1/(se.outcome^2 + bvb + tau2))
        1/p * t(beta.outcome - beta.exposure %*% beta.hat) %*% S %*%
          (beta.outcome - beta.exposure %*% beta.hat)})
  phi_selected <- phi_cand[which.min(prof.lik)]
  MV.l.inv <- MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_selected/MV_eigvalues)) %*% t(MV_eigen$vectors)
  mvmr.adIVW <- MV.l.inv %*% t(beta.exposure) %*% W %*% beta.outcome
  if (over.dispersion) {
    tau2_adivw <- Reduce("+",lapply(1:p, function(j) {
      v <- Vj[[j]] * (se.outcome[j]^(-2))
      (beta.outcome[j] - beta.exposure[j,] %*% mvmr.adIVW)^2*se.outcome[j]^(-2) - 1 - as.numeric(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
      }))/sum(diag(W))
    tau2_adivw <- as.numeric(tau2_adivw)
  } else {
    tau2_adivw <- 0
  }
  adIVW_Vt <- Reduce("+",lapply(1:p, function(j) {
    m <- beta.exposure[j,] %*% t(beta.exposure[j,]) * (se.outcome[j]^(-2))
    v <- Vj[[j]]*(se.outcome[j]^(-2))
    bvb <- as.numeric(t(mvmr.adIVW) %*% v %*% mvmr.adIVW)
    vbbv <- v %*% mvmr.adIVW %*% t(mvmr.adIVW) %*% v
    m*(1+bvb+tau2_adivw*se.outcome[j]^(-2)) + vbbv
  }))
  mvmr.adIVW.se <- sqrt(diag(MV.l.inv%*%adIVW_Vt%*%MV.l.inv))
  return(list(beta.hat = mvmr.adIVW,
              beta.se = mvmr.adIVW.se,
              iv_strength_parameter = iv_strength_parameter,
              phi_selected = phi_selected,
              tau.square = tau2_adivw))
}
