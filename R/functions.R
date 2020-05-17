Phi.tilde<-function(x){
  res<-pnorm(x,lower.tail = FALSE)
  return(res)
}

ivw<-function(beta.exposure, beta.outcome, se.exposure, se.outcome, alpha=0.05, pval.selection=NULL,lambda=0){
  if(lambda==0){
    ind<-1:length(beta.exposure)
  }else if (lambda>0){
    if(is.null(pval.selection)){stop("When lambda>0, need to provide pval.selection.")}
    else{
      pvalue<-2*Phi.tilde(lambda)
      ind<-which(pval.selection<=pvalue)
    }
  }
  beta_IVW<-sum(beta.outcome[ind]*beta.exposure[ind]/(se.outcome[ind])^2)/sum(beta.exposure[ind]^2/(se.outcome[ind])^2)
  se.ratio<-se.exposure/se.outcome
  mu<-beta.exposure/se.exposure
  V1<-sum((se.ratio^2*mu^2+beta_IVW^2*se.ratio^4*(mu^2+1))[ind])
  V2<-sum((se.ratio^2*mu^2)[ind])
  var_IVW<-V1/V2^2
  se_IVW<-sqrt(var_IVW)
  c_alpha<-qnorm(alpha/2,lower.tail = FALSE)
  CI_IVW.lower<-beta_IVW-c_alpha*se_IVW
  CI_IVW.upper<-beta_IVW+c_alpha*se_IVW
  return(list(beta.hat=beta_IVW,beta.se=se_IVW,n.IV=length(ind),IV=ind))
}

#' Main function for dIVW
#'
#' @param beta.exposure A vector of SNP effects on the exposure vairable, usually obtained from a GWAS
#' @param beta.outcome A vector of SNP effects on the outcome vairable, usually obtained from a GWAS
#' @param se.exposure A vecor of standard errors of \code{beta.exposure}
#' @param se.outcome A vector of standard errors of \code{beta.outcome}
#' @param alpha Confidence interval has level 1-alpha
#' @param pval.selection A vector of p-values calculated based on the selection dataset that is used for IV selection. It is not required when lambda=0
#' @param lambda The specified z-score threhold. Default is 0 (without thresholding)
#' @param over.dispersion Should the model consider balanced horizontal pleiotropy. Default is FALSE
#' @param diagnostics Should the function returns the q-q plot for assumption diagnosis. Default is FALSE
#'
#' @return A list
#' \describe{
#' \item{beta.hat}{Estimated causal effect}
#' \item{beta.se}{Standard error of \code{beta.hat}}
#' \item{condition}{A measure that needs to be large for reliable asymptotic approximation based on the dIVW estimator. It is recommended to be greater than 20}
#' \item{tau.square}{Overdispersion parameter if \code{over.dispersion=TRUE}}
#' \item{n.IV}{Number of IVs used in the dIVW estimator}
#' \item{IV}{IVs that are used in the dIVW estimator}
#' }
#'
#' @references Ting Ye, Jun Shao, Hyunseung Kang (2020). Debiased Inverse-Variance Weighted Estimator in Two-Sample Summary-Data Mendelian Randomization.\url{https://arxiv.org/abs/1911.09802}.
#'
#' @import stats
#' @export
#'
#' @examples
#'
#' data(bmi.cad)
#' attach(bmi.cad)
#' mr.divw(beta.exposure, beta.outcome, se.exposure, se.outcome, diagnostics=TRUE)
#' detach(bmi.cad)
#'
mr.divw<-function(beta.exposure, beta.outcome, se.exposure, se.outcome, alpha=0.05, pval.selection=NULL,lambda=0, over.dispersion=FALSE,diagnostics=FALSE){
  if(lambda==0){
    ind<-1:length(beta.exposure)
  }else if (lambda>0){
    if(is.null(pval.selection)){stop("When lambda>0, need to provide pval.selection.")}
    else{
      pvalue<-2*Phi.tilde(lambda)
      ind<-which(pval.selection<=pvalue)
    }
  }
  beta_dIVW<-sum(beta.outcome[ind]*beta.exposure[ind]/(se.outcome[ind])^2)/sum((beta.exposure[ind]^2-se.exposure[ind]^2)/(se.outcome[ind])^2)
  se.ratio<-se.exposure/se.outcome
  mu<-beta.exposure/se.exposure
  condition<-(mean((mu[ind])^2)-1)*sqrt(length(ind))/max(1,lambda^2)
  tau.square<-0
  if(over.dispersion){
    beta_dIVW_0<-sum(beta.outcome*beta.exposure/(se.outcome)^2)/sum((beta.exposure^2-se.exposure^2)/(se.outcome)^2)
    tau.square<-sum(((beta.outcome-beta_dIVW_0*beta.exposure)^2-se.outcome^2-beta_dIVW_0^2*se.exposure^2)/se.outcome^2)/sum(se.outcome^(-2))
  }
  V1<-sum((se.ratio^2*mu^2+beta_dIVW^2*se.ratio^4*(mu^2+1)+tau.square*se.ratio^2/se.outcome^2*mu^2)[ind])
  V2<-sum((se.ratio^2*(mu^2-1))[ind])
  var_dIVW<-V1/V2^2
  se_dIVW<-sqrt(var_dIVW)
  c_alpha<-qnorm(alpha/2,lower.tail = FALSE)
  CI_dIVW.lower<-beta_dIVW-c_alpha*se_dIVW
  CI_dIVW.upper<-beta_dIVW+c_alpha*se_dIVW
  if(diagnostics){
    t<-(beta.outcome-beta_dIVW*beta.exposure)/sqrt(se.outcome^2+tau.square+beta_dIVW^2*se.exposure^2)
    qqnorm(t)
    qqline(t)
  }
  return(list(beta.hat=beta_dIVW,beta.se=se_dIVW,condition=condition,tau.square=tau.square,n.IV=length(ind),IV=ind))
}

var.divw<-function(lambda,pval.selection, beta, se.ratio, mu, tau.square, se.outcome){
  pvalue<-2*Phi.tilde(lambda)
  ind<-which(pval.selection<=pvalue)
  V1<-sum((se.ratio^2*mu^2+beta^2*se.ratio^4*(mu^2+1)+tau.square*se.ratio^2/se.outcome^2*mu^2)[ind])
  V2<-sum((se.ratio^2*(mu^2-1))[ind])
  res<-V1/V2^2
  return(res)
}

#' MR-EO Algorithm to Adaptively Find the Optimal Z-score Threhold.
#'
#' @param lambda.start Initial value for lambda (the z-score threshold).
#' @param beta.exposure A vector of SNP effects on the exposure vairable, usually obtained from a GWAS
#' @param beta.outcome A vector of SNP effects on the outcome vairable, usually obtained from a GWAS
#' @param se.exposure A vecor of standard errors of \code{beta.exposure}
#' @param se.outcome A vector of standard errors of \code{beta.outcome}
#' @param pval.selection A vector of p-values calculated based on the selection dataset that is used for IV selection
#' @param over.dispersion Should the model consider balanced horizontal pleiotropy. Default is FALSE
#' @param max_opt_iter Maximum number of iterations. Default is 5
#'
#' @details \code{mr.eo} is an adaptive algorithm that finds the optimal z-socre threshold that leads to the dIVW estimator with the smallest variance.
#' @return A list
#' \describe{
#' \item{lambda.opt}{Optimal z-socre threshold}
#' \item{n.iter}{Number of iterations to find \code{lambda.opt}}
#' }
#'
#' @references Ting Ye, Jun Shao, Hyunseung Kang (2020). Debiased Inverse-Variance Weighted Estimator in Two-Sample Summary-Data Mendelian Randomization.\url{https://arxiv.org/abs/1911.09802}.
#'
#' @import stats
#' @export
#' @examples
#'
#' df<-data_gen_summary("case1")
#' attach(df)
#' lambda.opt<-mr.eo(0, beta.exposure, beta.outcome, se.exposure, se.outcome, pval.selection)$lambda.opt
#' mr.divw(beta.exposure, beta.outcome, se.exposure, se.outcome, pval.selection=pval.selection, lambda=lambda.opt)
#' detach(df)
#'
#' data(bmi.cad)
#' attach(bmi.cad)
#' lambda.opt<-mr.eo(0, beta.exposure, beta.outcome, se.exposure, se.outcome, pval.selection)$lambda.opt
#' mr.divw(beta.exposure, beta.outcome, se.exposure, se.outcome, pval.selection=pval.selection, lambda=lambda.opt)
#' detach(bmi.cad)
#'
mr.eo<-function(lambda.start, beta.exposure, beta.outcome, se.exposure, se.outcome, pval.selection, over.dispersion=FALSE, max_opt_iter=5){
  beta_res<-numeric(max_opt_iter)
  lambda_res<-numeric(max_opt_iter)
  var_res<-numeric(max_opt_iter)
  se.ratio<-se.exposure/se.outcome
  mu<-beta.exposure/se.exposure
  max_lambda<-qnorm(min(pval.selection)/2,lower.tail = FALSE)
  opt_right_end<-min(max_lambda,sqrt(2*log(length(beta.exposure))))
  iter<-1
  lambda_res[iter]<-lambda.start
  if(lambda.start>=max_lambda){
    warning("The initial lambda is too large and fails to select any IV. The initial lambda is reset to 0.")
    lambda_res[iter]<-0
  }
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection=pval.selection,
                    lambda=lambda_res[iter], over.dispersion = over.dispersion)
  beta_res[iter]<-tmp$beta.hat
  tau.square<-tmp$tau.square
  var_res[iter]<-var.divw(lambda_res[iter],pval.selection,beta_res[iter],se.ratio,mu,tau.square,se.outcome)
  lambda_res[iter+1]<-optimize(var.divw,c(0,opt_right_end),tol=0.001,pval.selection=pval.selection,beta=beta_res[iter],
                               se.ratio=se.ratio,mu=mu, tau.square=tau.square,se.outcome=se.outcome)$minimum
  iter<-2
  beta_res[iter]<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection=pval.selection,
                       lambda=lambda_res[iter], over.dispersion = over.dispersion)$beta.hat
  var_res[iter]<-var.divw(lambda_res[iter],pval.selection,beta_res[iter],se.ratio,mu,tau.square,se.outcome)
  while(iter<=max_opt_iter & var_res[iter]<var_res[iter-1]){
    lambda_res[iter+1]<-optimize(var.divw,c(0,opt_right_end),tol=0.001,pval.selection=pval.selection,beta=beta_res[iter],
                                 se.ratio=se.ratio,mu=mu, tau.square=tau.square,se.outcome=se.outcome)$minimum
    iter<-iter+1
    beta_res[iter]<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection=pval.selection,
                         lambda=lambda_res[iter], over.dispersion = over.dispersion)$beta.hat
    var_res[iter]<-var.divw(lambda_res[iter],pval.selection,beta_res[iter],se.ratio,mu,tau.square,se.outcome)
  }
  return(list(lambda.opt=lambda_res[iter-1],n.iter=iter-1))
}

#' Summary Statistics Simulated from the BMI-CAD Dataset
#'
#' @param case Simulation scenario used in Section 5.1 in Ye et al., (2020).
#'
#' @return A data frame
#' @references Ting Ye, Jun Shao, Hyunseung Kang (2020). Debiased Inverse-Variance Weighted Estimator in Two-Sample Summary-Data Mendelian Randomization.\url{https://arxiv.org/abs/1911.09802}.
#'
#' @import stats
#' @export
#'
data_gen_summary<-function(case=c("case1","case2","case3","case3_pleiotropy")){
  p<-1119
  beta0<-0.4
  data(bmi.cad)
  pi<-bmi.cad$beta.exposure
  if(case=="case1"){
    s<-20
    strong_ind<-order((bmi.cad$beta.exposure/bmi.cad$se.exposure)^2,decreasing = TRUE)[1:s]
    strong_pi<-pi[strong_ind]
    set.seed(0) # this seed is for null ses
    null_se.exposure<-sample(bmi.cad$se.exposure,replace = TRUE,size = p-s)
    null_se.outcome<-sample(bmi.cad$se.outcome,replace = TRUE,size=p-s)
    null_se.selection<-sample(bmi.cad$se.selection,replace = TRUE,size=p-s)
    useful_df<-data.frame(se.exposure=bmi.cad$se.exposure[strong_ind],
                          se.outcome=bmi.cad$se.outcome[strong_ind],
                          se.selection=bmi.cad$se.selection[strong_ind])
    useful_df$beta.exposure<-rnorm(s,mean=strong_pi,sd = useful_df$se.exposure)
    useful_df$beta.outcome<-rnorm(s,mean=beta0*strong_pi,sd=useful_df$se.outcome)
    useful_df$beta.selection<-rnorm(s,mean=strong_pi,sd = useful_df$se.selection)
    useful_df$relevant.ind<-1
    null_df<-data.frame(se.exposure=null_se.exposure,se.outcome=null_se.outcome,se.selection=null_se.selection)
    null_df$beta.exposure<-rnorm(p-s,0,sd=null_df$se.exposure)
    null_df$beta.outcome<-rnorm(p-s,0,sd=null_df$se.outcome)
    null_df$beta.selection<-rnorm(p-s,0,sd=null_df$se.selection)
    null_df$relevant.ind<-0
    full_df<-rbind(useful_df,null_df)
    full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
    full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
    full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
  }else if (case=="case2"){
    s<-100
    weak_ind<-1:s
    weak_pi<-pi[weak_ind]
    set.seed(0) # this seed is for null ses
    null_se.exposure<-sample(bmi.cad$se.exposure,replace = TRUE,size = p-s)
    null_se.outcome<-sample(bmi.cad$se.outcome,replace = TRUE,size=p-s)
    null_se.selection<-sample(bmi.cad$se.selection,replace = TRUE,size=p-s)
    useful_df<-data.frame(se.exposure=bmi.cad$se.exposure[weak_ind],
                          se.outcome=bmi.cad$se.outcome[weak_ind],
                          se.selection=bmi.cad$se.selection[weak_ind])
    useful_df$beta.exposure<-rnorm(s,mean=weak_pi,sd = useful_df$se.exposure)
    useful_df$beta.outcome<-rnorm(s,mean=beta0*weak_pi,sd=useful_df$se.outcome)
    useful_df$beta.selection<-rnorm(s,mean=weak_pi,sd = useful_df$se.selection)
    useful_df$relevant.ind<-1
    null_df<-data.frame(se.exposure=null_se.exposure,se.outcome=null_se.outcome,se.selection=null_se.selection)
    null_df$beta.exposure<-rnorm(p-s,0,sd=null_df$se.exposure)
    null_df$beta.outcome<-rnorm(p-s,0,sd=null_df$se.outcome)
    null_df$beta.selection<-rnorm(p-s,0,sd=null_df$se.selection)
    null_df$relevant.ind<-0
    full_df<-rbind(useful_df,null_df)
    full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
    full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
    full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
  }else if (case %in% c("case3","case3_pleiotropy")){
    full_df<-data.frame(se.exposure=bmi.cad$se.exposure,
                        se.outcome=bmi.cad$se.outcome,
                        se.selection=bmi.cad$se.selection)
    full_df$beta.exposure<-rnorm(p,mean=pi,sd=full_df$se.exposure)
    if(case=="case3_pleiotropy"){
      alpha<-rnorm(p,0,2/p*sum(full_df$se.outcome))
    }else{
      alpha<-rep(0,p)
    }
    full_df$beta.outcome<-rnorm(p,mean=pi*beta0+alpha,sd=full_df$se.outcome)
    full_df$beta.selection<-rnorm(p,mean=pi,sd=full_df$se.selection)
    full_df$relevant.ind<-1
    full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
    full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
    full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
  }
  return(full_df)
}

#' Summary Statistics from Simulated Individual-level data
#'
#' @param case Simulation scenario used in Section 5.3 in Ye et al., (2020).
#' @param true_var Do se.exposure, se.outcome, se.selection equal the true standard deviations or the estimated standard errors. Default is FALSE.
#'
#' @return A data frame
#' @references Ting Ye, Jun Shao, Hyunseung Kang (2020). Debiased Inverse-Variance Weighted Estimator in Two-Sample Summary-Data Mendelian Randomization.\url{https://arxiv.org/abs/1911.09802}.
#'
#' @import stats
#' @export
#'
data_gen_individual<-function(case=c("case4","case5","case6","case7"),true_var=FALSE){
  if(case=="case4"){
    s<-200
    n<-1e4
  }else if (case=="case5"){
    s<-1e3
    n<-1e4
  }else if (case=="case6"){
    s<-1e3
    n<-5e4
  }else if (case=="case7"){
    s<-2e3
    n<-1e4
  }
  p<-2000
  h<-sqrt(0.2) # 20% total heritability
  beta0<-0.4
  standard_normal_random_effect<-mr.divw:::standard_normal_random_effect
  gamma_coef<-standard_normal_random_effect*(h/sqrt(s))*sqrt(2)
  if(s<p){
    gamma_coef[(s+1):p]<-0
  }
  var_x<-0.5*sum(gamma_coef[1:s]^2 )+(1-h^2)
  var_y<-beta0^2*var_x+(2*beta0+1)*(1-h^2)*0.6+1
  sigma_yj_true<-sqrt((n*0.5)^(-1)*(var_y-beta0^2*gamma_coef[1:p]^2*0.5))
  sigma_xj_true<-sqrt((n*0.5)^(-1)*(var_x-gamma_coef[1:p]^2*0.5))
  data_gen_onesample<-function(){
    z<-matrix(nrow=n,ncol=p)
    for(j in 1:p){
      tmp<-rmultinom(n,1,c(0.25,0.5,0.25)) # z's are independent
      z[,j]<-0*(tmp[1,]==1)+1*(tmp[2,]==1)+2*(tmp[3,]==1)
    }
    u<-rnorm(n,0,sqrt((1-h^2)*0.6))
    x<-z[,1:s] %*% gamma_coef[1:s] +u+rnorm(n,0,sqrt((1-h^2)*0.4))
    y<-10+beta0*x+u+rnorm(n)
    return(list(z=z,x=x,y=y))
  }
  tmp1<-data_gen_onesample() # exposure dataset
  tmp2<-data_gen_onesample() # outcome dataset
  tmp3<-data_gen_onesample() # selection dataset
  full_df<-data.frame(beta.exposure=numeric(p))
  full_df$beta.exposure<-sapply(1:p, function(i) {
    coef <- lm(tmp1$x ~ tmp1$z[,i])$coef[-1]
    return(coef)
  })
  full_df$beta.outcome<-sapply(1:p, function(i) {
    coef <- lm(tmp2$y ~ tmp2$z[,i])$coef[-1]
    return(coef)
  })
  full_df$beta.selection<-sapply(1:p, function(i) {
    coef <- lm(tmp3$x ~ tmp3$z[,i])$coef[-1]
    return(coef)
  })
  full_df$relevant.ind<-0
  full_df$relevant.ind[1:s]<-1
  #####################
  if(true_var==TRUE){
    full_df$se.exposure<-sigma_xj_true # true sd
    full_df$se.outcome<-sigma_yj_true  # true sd
    full_df$se.selection<-sigma_xj_true # true sd
  }else{
    full_df$se.exposure<-sapply(1:p, function(i) {
      se <- summary( lm(tmp1$x ~ tmp1$z[,i]))$coefficients[2,2]
      return(se)
    })
    full_df$se.outcome<-sapply(1:p, function(i) {
      se <- summary( lm(tmp2$y ~ tmp2$z[,i]))$coefficients[2,2]
      return(se)
    })
    full_df$se.selection<-sapply(1:p, function(i) {
      se <- summary(lm(tmp3$x ~ tmp3$z[,i]))$coefficients[2,2]
      return(se)
    })
  }
  #####################
  full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
  full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
  full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
  return(full_df)
}

#' Published Table 4
#'
#' @return A table
#'
#' @export
#' @examples
#' table_publish()
#'
table_publish<-function(){
  res<-matrix(data=NA, nrow = 6, ncol=5)
  rownames(res)<-c("lambda","n_IV","condition","IVW","dIVW","dIVW_alpha")
  res[1,1:3]<-c(0, 5.45, round(sqrt(2*log(1119)),2))
  data("bmi.cad")
  attach(bmi.cad)
  tmp<-ivw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=0)
  res[4,1]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[2,1]<-tmp$n.IV
  tmp<-ivw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=5.45)
  res[4,2]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[2,2]<-tmp$n.IV
  tmp<-ivw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=sqrt(2*log(1119)))
  res[4,3]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[2,3]<-tmp$n.IV
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=0)
  res[5,1]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[3,1]<-round(tmp$condition,1)
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=5.45)
  res[5,2]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[3,2]<-round(tmp$condition,1)
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=sqrt(2*log(1119)))
  res[5,3]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[3,3]<-round(tmp$condition,1)
  lambda.opt<-mr.eo(sqrt(2*log(1119)),beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection)$lambda.opt
  res[1,4]<-round(lambda.opt,2)
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=lambda.opt)
  res[5,4]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[3,4]<-round(tmp$condition,1)
  res[2,4]<-tmp$n.IV
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=0,over.dispersion = TRUE)
  res[6,1]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=5.45,over.dispersion = TRUE)
  res[6,2]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=sqrt(2*log(1119)),over.dispersion = TRUE)
  res[6,3]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  lambda.opt<-mr.eo(sqrt(2*log(1119)),beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection,TRUE)$lambda.opt
  res[1,5]<-round(lambda.opt,2)
  tmp<-mr.divw(beta.exposure,beta.outcome,se.exposure,se.outcome,pval.selection = pval.selection,lambda=lambda.opt,over.dispersion = TRUE)
  res[6,4]<-paste0(round(tmp$beta.hat,3)," (",round(tmp$beta.se,3),")")
  res[2,5]<-tmp$n.IV
  res[3,5]<-round(tmp$condition,1)
  detach(bmi.cad)
  return(res)
}

