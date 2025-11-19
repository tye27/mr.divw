#'  Perform spectral regularized inverse-variance weighted (SRIVW) estimator for summary-data multivariable Mendelian randomization
#'
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param phi_cand A vector of tuning parameters for SRIVW estimator. Default is 0. To use the recommended set for the tuning parameter, simply set phi_cand = NULL.
#' @param over.dispersion Should the model consider balanced horizontal pleiotropy? Default is TRUE.
#' @param overlap Should the model consider overlapping exposure and outcome datasets? Default is FALSE.
#' @param gen_cor If overlap = FALSE, provide a K-by-K matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure. If overlap = TRUE, provide a (K+1)-by-(K+1) matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure and the outcome, where the last index position corresponds to the outcome.  The correlations can either be estimated, be assumed to be zero, or fixed at zero. Default input is NULL, meaning that an identity matrix is used as the correlation matrix.
#'
#' @return A list with elements
#' \item{beta.hat}{Estimated direct effects of each exposure on the outcome}
#' \item{beta.se}{Estimated standard errors of beta.hat}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies the IV strength in the sample}
#' \item{phi_selected}{The selected tuning parameter for the SRIVW estimator}
#' \item{tau.square}{Overdispersion parameter if \code{over.dispersion=TRUE}}
#' @export
#'
#' @examples
#' data("hdl_subfractions")
#' # We are going to estimate the effect of S-HDL-P on the risk of CAS with adjustments of HDL, LDL, and TG levels
#' beta.exposure <- hdl_subfractions$data[,c("gamma_exp1","gamma_exp2","gamma_exp3","gamma_exp4")] # corresponding to SNP effects on respectively  HDL, LDL, TG, and S-HDL-P
#' se.exposure <- hdl_subfractions$data[,c("se_exp1","se_exp2","se_exp3","se_exp4")]
#' beta.outcome <- hdl_subfractions$data$gamma_out1
#' se.outcome <- hdl_subfractions$data$se_out1
#' P <- hdl_subfractions$cor.mat[c(1:4,10),c(1:4,10)] # make sure the last index corresponds to the outcome
#' mvmr.srivw(beta.exposure = beta.exposure,
#' se.exposure = se.exposure,
#' beta.outcome = beta.outcome,
#' se.outcome = se.outcome,
#' gen_cor = P,
#' phi_cand = NULL,
#' over.dispersion = FALSE,
#' overlap = TRUE)
#'
mvmr.srivw <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, phi_cand=0, over.dispersion = TRUE, overlap = FALSE, gen_cor = NULL) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {stop("this function is developed for multivariable MR")}
  K <- ncol(beta.exposure)
  if (is.null(gen_cor)) {
    if (overlap == FALSE) {P <- diag(K)} else {P <- diag(K+1)}
  } else {
    P <- as.matrix(gen_cor)
  }
  if (overlap == FALSE & (ncol(P) != ncol(beta.exposure))) {stop("The shared correlation matrix has a different number of columns than the input beta.exposure")}
  if (overlap == TRUE & (ncol(P) != (ncol(beta.exposure) + 1))) {stop("If the exposure and outcome datasets are overlapping, please provide a (K+1)-by-(K+1) correlation matrix.")}
  if (nrow(beta.exposure) != length(beta.outcome)) {stop("The number of SNPs in beta.exposure and beta.outcome is different")}
  if (over.dispersion == TRUE & overlap == TRUE) {stop("The function allowing for both balanced horizontal pleiotropy and overlapping datasets is still under development")}
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  # the number of SNPs
  p <- nrow(beta.exposure)
  # diagonal W matrix
  W<- diag(se.outcome^(-2))
  # create a list of Sigma Xj matrices
  if (overlap == FALSE) {
    Vj <- lapply(1:p, function(j) diag(se.exposure[j,]) %*% P %*% diag(se.exposure[j,]))
  } else {
    Vj <- lapply(1:p, function(j) diag(c(se.exposure[j,],se.outcome[j])) %*% P %*% diag(c(se.exposure[j,],se.outcome[j])))
  }
  # calculate square root inverse of P
  P_eigen <- eigen(P[1:K,1:K])
  P_root_inv <- P_eigen$vectors %*% diag(1/sqrt(P_eigen$values)) %*% t(P_eigen$vectors)
  Vj_root_inv <- lapply(1:p, function(j) {
    P_root_inv %*% diag(1/se.exposure[j,])
  })
  # calcualte IV strenght parameter
  IV_strength_matrix <- Reduce("+",lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j,]
    beta.exposure.V %*% t(beta.exposure.V)})) - p*diag(K)
  iv_strength_parameter <- min(eigen(IV_strength_matrix/sqrt(p))$values)
  if (iv_strength_parameter < 7) {
    warning("Sample IV strength parameter is too small, indicating potential weak instrument bias and incorrect inference.")
  }
  # get V matrix
  V <- Reduce("+",lapply(1:p, function(j) {Vj[[j]][1:K,1:K] * (se.outcome[j]^(-2))}))
  # get M matrix
  M <- t(beta.exposure)%*%W%*%beta.exposure
  # get M-V matrix
  MV <- M-V
  # perform eigen-decomposition on MV
  MV_eigvalues <- eigen(MV)$values
  MV_eigen <- eigen(MV)
  # multivariable (a)dIVW estimator
  if (is.null(phi_cand)) {
    phi_cand <- c(0, exp(seq(0, 17, by = 0.5) - min(iv_strength_parameter)))
  }
  phi_length <- length(phi_cand)
  MV.l.inv.long <- Reduce(rbind, lapply(1:phi_length, function(l) {
    MV_eigen$vectors %*% diag(1/(MV_eigvalues + phi_cand[l]/MV_eigvalues)) %*% t(MV_eigen$vectors)}
  ))
  if (overlap == FALSE) {# (a)dIVW for independent datasets
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
      if (tau2_adivw < 0) {
        tau2_adivw <- 0
        warning("Estimated overdispersion parameter < 0. Fixed at 0 instead.")
      }
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
  } else {
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
    tau2_adivw <- NULL # current version does not allow for over-dispersion and overlapping datasets at the same time.
  }
  return(list(beta.hat = mvmr.adIVW,
              beta.se = mvmr.adIVW.se,
              iv_strength_parameter = iv_strength_parameter,
              phi_selected = phi_selected,
              tau.square = tau2_adivw))
}
