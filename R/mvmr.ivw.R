#' Perform inverse-variance weighted (IVW) estimator for two-sample summary-data multivariable Mendelian randomization
#'
#' @param beta.exposure A data.frame or matrix. Each row contains the estimated marginal effect of a SNP on K exposures, usually obtained from a GWAS
#' @param se.exposure A data.frame or matrix of estimated standard errors of beta.exposure
#' @param beta.outcome A vector of the estimated marginal effect of a SNP on outcome, usually obtained from a GWAS
#' @param se.outcome A vector of estimated standard errors of beta.outcome
#' @param gen_cor A K-by-K matrix for the estimated shared correlation matrix between the effect of the genetic variants on each exposure, where K is the number of exposure. The correlations can either be estimated, be assumed to be zero, or fixed at zero using non-overlapping samples of each exposure GWAS. Default input is NULL, meaning that an identity matrix is used as the correlation matrix.
#'
#' @return A list with elements
#' \item{beta.hat}{Estimated direct effects of each exposure on the outcome}
#' \item{beta.se}{Estimated standard errors of beta.hat}
#' \item{iv_strength_parameter}{The minimum eigenvalue of the sample IV strength matrix, which quantifies the IV strength in the sample}
#' @import MVMR
#' @export
#'
#' @examples
#' data("hdl_subfractions")
#' # We are going to estimate the effect of S-HDL-P on the risk of CAS with adjustments of HDL, LDL, and TG levels
#' beta.exposure <- hdl_subfractions$data[,c("gamma_exp1","gamma_exp2","gamma_exp3","gamma_exp4")]
#' se.exposure <- hdl_subfractions$data[,c("se_exp1","se_exp2","se_exp3","se_exp4")]
#' beta.outcome <- hdl_subfractions$data$gamma_out1
#' se.outcome <- hdl_subfractions$data$se_out1
#' P <- hdl_subfractions$cor.mat[c(1:4),c(1:4)]
#' mvmr.ivw(beta.exposure = beta.exposure,
#' se.exposure = se.exposure,
#' beta.outcome = beta.outcome,
#' se.outcome = se.outcome,
#' gen_cor = P)
#'
mvmr.ivw <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = NULL) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {stop("either beta.exposure or se.exposure only has one column; if univariable MR is performed, please use univariable MR methods.")}
  # the number of exposures
  K <- ncol(beta.exposure)
  if (is.null(gen_cor)) {P <- diag(K)} else {P <- as.matrix(gen_cor)}
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
  # multivariable IVW estimator
  mvmr.IVW <- solve(M) %*% t(beta.exposure) %*% W %*% beta.outcome
  # calculate variance covariance matrix
  IVW_Vt <- Reduce("+", lapply(1:p, function(j) {
    m <- beta.exposure[j,] %*% t(beta.exposure[j,]) * (se.outcome[j]^(-2))
    v <- Vj[[j]] * (se.outcome[j]^(-2))
    bvb <- as.numeric(t(mvmr.IVW) %*% v %*% mvmr.IVW)
    vbbv <- v %*% mvmr.IVW %*% t(mvmr.IVW) %*% v
    m*(1+bvb) + vbbv
  }))
  mvmr.IVW.se <- sqrt(diag(solve(M)%*%IVW_Vt%*%solve(M)))
  return(list(beta.hat = mvmr.IVW,
              beta.se = mvmr.IVW.se,
              iv_strength_parameter = iv_strength_parameter))
}
