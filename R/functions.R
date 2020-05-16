Phi.tilde<-function(x){
  res<-pnorm(x,lower.tail = FALSE)
  return(res)
}

IVW<-function(df,lambda){
  pvalue<-2*Phi.tilde(lambda)
  ind<-which(df$pval.selection<=pvalue)
  beta_IVW<-sum(df$beta.outcome[ind]*df$beta.exposure[ind]/(df$se.outcome[ind])^2)/sum(df$beta.exposure[ind]^2/(df$se.outcome[ind])^2)
  df$se.ratio<-df$se.exposure/df$se.outcome
  df$mu<-df$beta.exposure/df$se.exposure
  condition<-(mean((df$mu[ind])^2)-1)*sqrt(length(ind))/max(1,lambda^2)
  var_IVW<-var_IVW_f(lambda,beta_IVW,df,var.method)
  se_IVW<-sqrt(var_IVW)
  c_alpha<-qnorm(0.025,lower.tail = FALSE)
  CI_IVW.lower<-beta_IVW-c_alpha*se_IVW
  CI_IVW.upper<-beta_IVW+(c_alpha)*se_IVW
  ind_CP<-1*(beta0>CI_IVW.lower & beta0 <CI_IVW.upper)
  return(list(n_IV=length(ind),IV=ind,beta=beta_IVW,se_beta=se_IVW,ind_CP=ind_CP, condition=condition))
}

var_IVW_f<-function(lambda,beta,df){
  pvalue<-2*Phi.tilde(lambda)
  ind<-which(df$pval.selection<=pvalue)
  df$se.ratio<-df$se.exposure/df$se.outcome
  df$mu<-df$beta.exposure/df$se.exposure
  V1<-sum((df$se.ratio^2*df$mu^2+beta^2*df$se.ratio^4*(df$mu^2+1))[ind])
  V2<-sum((df$se.ratio^2*df$mu^2)[ind])
  var_IVW<-V1/V2^2
  return(var_IVW)
}

dIVW<-function(df,lambda=0,over.dispersion=FALSE,var.method=c("bootstrap","empirical")){
  pvalue<-2*Phi.tilde(lambda)
  ind<-which(df$pval.selection<=pvalue)
  beta_dIVW<-sum(df$beta.outcome[ind]*df$beta.exposure[ind]/(df$se.outcome[ind])^2)/sum((df$beta.exposure[ind]^2-df$se.exposure[ind]^2)/(df$se.outcome[ind])^2)
  df$se.ratio<-df$se.exposure/df$se.outcome
  df$mu<-df$beta.exposure/df$se.exposure
  condition<-(mean((df$mu[ind])^2)-1)*sqrt(length(ind))/max(1,lambda^2)
  var_dIVW<-var_dIVW_f(lambda,beta_dIVW,df,over.dispersion,var.method)
  se_dIVW<-sqrt(var_dIVW)
  c_alpha<-qnorm(0.025,lower.tail = FALSE)
  CI_dIVW.lower<-beta_dIVW-c_alpha*se_dIVW
  CI_dIVW.upper<-beta_dIVW+(c_alpha)*se_dIVW
  CI_len<-2*c_alpha*se_dIVW
  ind_CP<-1*(beta0>CI_dIVW.lower & beta0 <CI_dIVW.upper)
  return(list(n_IV=length(ind),IV=ind,beta=beta_dIVW,se_beta=se_dIVW, CI_len=CI_len,ind_CP=ind_CP, condition=condition))
}

var_dIVW_f<-function(lambda,beta,df,over.dispersion){
  pvalue<-2*Phi.tilde(lambda)
  ind<-which(df$pval.selection<=pvalue)
  df$se.ratio<-df$se.exposure/df$se.outcome
  df$mu<-df$beta.exposure/df$se.exposure
  tau.square<-0
  if(over.dispersion==TRUE){
      beta_dIVW_0<-sum(df$beta.outcome*df$beta.exposure/(df$se.outcome)^2)/sum((df$beta.exposure^2-df$se.exposure^2)/(df$se.outcome)^2)
      tau.square<-sum(((df$beta.outcome-beta_dIVW_0*df$beta.exposure)^2-df$se.outcome^2-beta_dIVW_0^2*df$se.exposure^2)/df$se.outcome^2)/sum(df$se.outcome^(-2))
  }
  V1<-sum((df$se.ratio^2*df$mu^2+beta^2*df$se.ratio^4*(df$mu^2+1)+tau.square*df$se.ratio^2/df$se.outcome^2*df$mu^2)[ind])
  V2<-sum((df$se.ratio^2*(df$mu^2-1))[ind])
  var_dIVW<-V1/V2^2
  return(var_dIVW)
}

opt.lambda.alternating<-function(df,lambda.start,over.dispersion,var.method=c("bootstrap","empirical"),max_opt_iter=10){
  opt_right_end<-min(max(df$beta.selection/df$se.selection),sqrt(2*log(dim(df)[1])))
  initial.beta<-dIVW(df,lambda.start, over.dispersion,"empirical")$beta_dIVW
  beta_res<-numeric(max_opt_iter)
  lambda_res<-numeric(max_opt_iter)
  var_res<-numeric(max_opt_iter)
  iter<-1
  lambda_res[iter]<-lambda.start
  tmp_n_IV<-dIVW(df,lambda_res[iter],over.dispersion,var.method)$n_IV
  if(tmp_n_IV==0){
    warning("The initial lambda is too large and fails to select any IV. The initial lambda is reset to 0.")
    lambda_res[iter]<-0
  }
  beta_res[iter]<-dIVW(df,lambda_res[iter],over.dispersion, var.method)$beta
  var_res[iter]<-var_dIVW_f(lambda_res[iter],beta_res[iter],df,over.dispersion,var.method)
  lambda_res[iter+1]<-optimize(var_dIVW_f,c(0,opt_right_end),tol=0.001,beta=beta_res[iter],df=df,over.dispersion=over.dispersion,var.method=var.method)$minimum
  iter<-2
  beta_res[iter]<-dIVW(df,lambda_res[iter],over.dispersion, var.method)$beta
  var_res[iter]<-var_dIVW_f(lambda_res[iter],beta_res[iter],df,over.dispersion,var.method)
  while(iter<=max_opt_iter & var_res[iter]<var_res[iter-1]){
    lambda_res[iter+1]<-optimize(var_dIVW_f,c(0,opt_right_end),tol=0.001,beta=beta_res[iter],df=df,over.dispersion=over.dispersion,var.method=var.method)$minimum
    iter<-iter+1
    beta_res[iter]<-dIVW(df,lambda_res[iter],over.dispersion, var.method)$beta
    var_res[iter]<-var_dIVW_f(lambda_res[iter],beta_res[iter],df,over.dispersion,var.method)
  }
  return(list(lambda.iter=lambda_res[iter-1],n.iter=iter-1))
}


diagostic_plot<-function(df,over.dispersion){
  tmp<-IVW_subset(df,over.dispersion = over.dispersion)
  t<-(df$beta.outcome-tmp$beta_IVW*df$beta.exposure)/sqrt(df$se.outcome^2+tmp$tau.square_est+tmp$beta_IVW^2*df$se.exposure^2)
  qqnorm(t)
  qqline(t)
}




