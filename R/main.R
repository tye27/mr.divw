##### Condor ######
##### Last update by Ting Ye, Mar 27, 2020 #####
#setwd("~/Dropbox (Penn)/!! My experience/Research/IV/Simulation/2020 Mar revision")
arg<-commandArgs(trailingOnly = TRUE)
njob<-eval(parse(text = arg[[1]]))
library(mr.raps)
#library(MendelianRandomization)
#library(rmutil)
source("functions.R")
source("setting1.R")

n_rep<-100
c_alpha<-qnorm(0.025,lower.tail = FALSE)
max_opt_iter<-5 
lambda_our<-sqrt(2*log(p))
lambda_5.45<-qnorm(2.5e-08,0,1,lower.tail = FALSE)
methods_name<-c("IVW_P","IVW_5.45","IVW_our","dIVW_P","dIVW_5.45","dIVW_our","dIVW_iter","IVW_P_boot",
                "IVW_5.45_boot","IVW_our_boot","dIVW_P_boot","dIVW_5.45_boot","dIVW_our_boot","dIVW_iter_boot",
                "raps","raps.shrink")
n_methods<-length(methods_name)

res<-matrix(nrow=n_rep,ncol=n_methods)
res_se<-matrix(nrow=n_rep,ncol=n_methods)
res_CP<-matrix(0,nrow=1,ncol=n_methods)
colnames(res)<-methods_name
colnames(res_se)<-methods_name
colnames(res_CP)<-methods_name
res_inclusion_our<-matrix(0,nrow=n_rep,ncol=(s+1))
res_inclusion_5.45<-matrix(0,nrow=n_rep,ncol=(s+1))
res_inclusion_opt<-matrix(0,nrow=n_rep,ncol=(s+1))
res_inclusion_opt_boot<-matrix(0,nrow=n_rep,ncol=(s+1))
res_lambda<-numeric(n_rep)
res_n_iter<-numeric(n_rep)
IV_names<-c(paste0("IV",1:s),"# unuseful IV")
colnames(res_inclusion_our)<-IV_names
colnames(res_inclusion_5.45)<-IV_names
colnames(res_inclusion_opt)<-IV_names
colnames(res_inclusion_opt_boot)<-IV_names


for(k in 1:n_rep){
  set.seed((njob-1)*n_rep+k)
  print(k)
  df_simu<-data_gen(true_var)
  
  var.method<-"empirical"
  #### IVW ####
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-IVW(df_simu,0,var.method = var.method) 
  res[k,1]<-tmp$beta
  res_se[k,1]<-tmp$se_beta
  res_CP[1]<-res_CP[1]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-IVW(df_simu,lambda_5.45,var.method=var.method)
  res[k,2]<-tmp$beta
  res_se[k,2]<-tmp$se_beta
  res_CP[2]<-res_CP[2]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-IVW(df_simu,lambda_our,var.method=var.method)
  res[k,3]<-tmp$beta
  res_se[k,3]<-tmp$se_beta
  res_CP[3]<-res_CP[3]+tmp$ind_CP
  
  #### dIVW ####
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=0, over.dispersion=over.dispersion, var.method=var.method)
  res[k,4]<-tmp$beta
  res_se[k,4]<-tmp$se_beta
  res_CP[4]<-res_CP[4]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=lambda_5.45,over.dispersion=over.dispersion, var.method=var.method)
  res[k,5]<-tmp$beta
  res_se[k,5]<-tmp$se_beta
  res_CP[5]<-res_CP[5]+tmp$ind_CP
  ind_5.45<-tmp$IV
  res_inclusion_5.45[k,1:s]<-1*((1:s)%in% ind_5.45)
  res_inclusion_5.45[k,s:(s+1)]<-length(which(!ind_5.45 %in% 1:s))
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=lambda_our,over.dispersion=over.dispersion, var.method=var.method)
  res[k,6]<-tmp$beta
  res_se[k,6]<-tmp$se_beta
  res_CP[6]<-res_CP[6]+tmp$ind_CP
  ind_our<-tmp$IV
  res_inclusion_our[k,1:s]<-1*((1:s)%in% ind_our)
  res_inclusion_our[k,s:(s+1)]<-length(which(!ind_our %in% 1:s))
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-opt.lambda.alternating(df = df_simu,lambda.start=sqrt(2*log(p)), 
                                   over.dispersion=over.dispersion,var.method=var.method,
                                   max_opt_iter = max_opt_iter)
  res_lambda[k]<-tmp$lambda.iter
  res_n_iter[k]<-tmp$n.iter
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=res_lambda[k],over.dispersion=over.dispersion, var.method=var.method)
  res[k,7]<-tmp$beta
  res_se[k,7]<-tmp$se_beta
  res_CP[7]<-res_CP[7]+tmp$ind_CP
  ind_opt<-tmp$IV
  res_inclusion_opt[k,1:s]<-1*((1:s)%in% ind_opt)
  res_inclusion_opt[k,s:(s+1)]<-length(which(!ind_opt %in% 1:s))
  
  var.method<-"bootstrap"
  #### IVW ####
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-IVW(df_simu,0,var.method = var.method) 
  res[k,8]<-tmp$beta
  res_se[k,8]<-tmp$se_beta
  res_CP[8]<-res_CP[8]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-IVW(df_simu,lambda_5.45,var.method=var.method)
  res[k,9]<-tmp$beta
  res_se[k,9]<-tmp$se_beta
  res_CP[9]<-res_CP[9]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-IVW(df_simu,lambda_our,var.method=var.method)
  res[k,10]<-tmp$beta
  res_se[k,10]<-tmp$se_beta
  res_CP[10]<-res_CP[10]+tmp$ind_CP
  
  #### dIVW ####
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=0,over.dispersion=over.dispersion, var.method=var.method)
  res[k,11]<-tmp$beta
  res_se[k,11]<-tmp$se_beta
  res_CP[11]<-res_CP[11]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=lambda_5.45,over.dispersion=over.dispersion, var.method=var.method)
  res[k,12]<-tmp$beta
  res_se[k,12]<-tmp$se_beta
  res_CP[12]<-res_CP[12]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=lambda_our,over.dispersion=over.dispersion, var.method=var.method)
  res[k,13]<-tmp$beta
  res_se[k,13]<-tmp$se_beta
  res_CP[13]<-res_CP[13]+tmp$ind_CP
  
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-opt.lambda.alternating(df = df_simu,lambda.start=sqrt(2*log(p)), 
                              over.dispersion=over.dispersion,var.method=var.method,
                              max_opt_iter = max_opt_iter)
  res_lambda[k]<-tmp$lambda.iter
  res_n_iter[k]<-tmp$n.iter
  tmp<-list(beta=NA,se_beta=NA,ind_CP=NA)
  tmp<-dIVW(df=df_simu,lambda=res_lambda[k],over.dispersion=over.dispersion, var.method=var.method)
  res[k,14]<-tmp$beta
  res_se[k,14]<-tmp$se_beta
  res_CP[14]<-res_CP[14]+tmp$ind_CP
  ind_opt<-tmp$IV
  res_inclusion_opt_boot[k,1:s]<-1*((1:s)%in% ind_opt)
  res_inclusion_opt_boot[k,s:(s+1)]<-length(which(!ind_opt %in% 1:s))
  
  
  
  # if(length(ind_our)<=2) {next}
  # df_our<-df_simu[ind_our,]
  # df_mr<-mr_input(bx = df_our$beta.exposure, bxse = df_our$se.exposure,
  #                 by = df_our$beta.outcome, byse = df_our$se.outcome)
  #  
  # 
  #  tmp10<-mr_egger(df_mr)
  #  res[k,10]<-tmp10$Estimate
  #  res_CI.len[k,10]<-2*c_alpha*tmp10$StdError.Est
  #  res_CP[10]<-res_CP[10]+1*(beta0>tmp10$CILower.Est &beta0< tmp10$CIUpper.Est)
  # 
  #  tmp11<-mr_median(df_mr)
  #  res[k,11]<-tmp11$Estimate
  #  res_CI.len[k,11]<-2*c_alpha*tmp11$StdError
  #  res_CP[11]<-res_CP[11]+1*(beta0>tmp11$CILower &beta0< tmp11$CIUpper)
  # 
  #  tmp12<-mr_mbe(df_mr)
  #  res[k,12]<-tmp12$Estimate
  #  res_CI.len[k,12]<-2*c_alpha*tmp12$StdError
  #  res_CP[12]<-res_CP[12]+1*(beta0>min(tmp12$CILower) &beta0< max(tmp12$CIUpper))
   
  
  #raps
  raps.fit1<-NA
  raps.fit1<-mr.raps(data=df_simu,over.dispersion = over.dispersion,loss.function = "l2",diagnostics = FALSE)
  if(!is.na(raps.fit1)){
    res[k,15]<-raps.fit1$beta.hat
    res_se[k,15]<-raps.fit1$beta.se
    res_CP[15]<-res_CP[15]+1*(beta0>(res[k,15]-c_alpha*raps.fit1$beta.se) &
                                beta0< (res[k,15]+c_alpha*raps.fit1$beta.se))
  }else{
    res[k,15]<-NA
    res_se[k,15]<-NA
    res_CP[15]<-NA
  }
  

  # #raps shrink
  raps.shrikage.fit<-NA
  z_tmp<-df_simu$beta.exposure/df_simu$se.exposure
  prior.par<-fit.mixture.model(z_tmp)
  raps.shrikage.fit<-mr.raps.shrinkage(b_exp=df_simu$beta.exposure,b_out = df_simu$beta.outcome,se_exp = df_simu$se.exposure,se_out = df_simu$se.outcome,loss.function = "l2",shrinkage = TRUE,prior.param = prior.par,over.dispersion = over.dispersion)
  if(!is.na(raps.shrikage.fit)){
    res[k,16]<-raps.shrikage.fit$beta.hat
    res_se[k,16]<-raps.shrikage.fit$beta.se
    res_CP[16]<-res_CP[16]+1*(beta0>(res[k,16]-c_alpha*raps.shrikage.fit$beta.se) &
                                beta0< (res[k,16]+c_alpha*raps.shrikage.fit$beta.se))
  }else{
    res[k,16]<-NA
    res_se[k,16]<-NA
    res_CP[16]<-NA
  }
  print(res[k,])
}

res<-as.data.frame(res)
write.csv(res,file=paste("Res/",file_name,njob,".csv",sep=""),row.names = FALSE)

res_se<-as.data.frame(res_se)
write.csv(res_se,file=paste("Res/",file_name,"_se_",njob,".csv",sep=""),row.names = FALSE)

res_CP<-as.data.frame(res_CP)
write.csv(res_CP,file=paste("Res/",file_name, "_CP_",njob,".csv",sep=""),row.names = FALSE)

write.csv(res_lambda, file=paste("Res/",file_name, "_lambda_",njob,".csv",sep=""), row.names = FALSE)
write.csv(res_n_iter, file=paste("Res/",file_name, "_n_iter_",njob,".csv",sep=""), row.names = FALSE)

res_inclusion_our<-as.data.frame(res_inclusion_our)
res_inclusion_5.45<-as.data.frame(res_inclusion_5.45)
res_inclusion_opt<-as.data.frame(res_inclusion_opt)
res_inclusion_opt_boot<-as.data.frame(res_inclusion_opt_boot)
write.csv(res_inclusion_our,file=paste("Res/",file_name,"_our_inclusion",njob,".csv",sep=""),row.names = FALSE)
write.csv(res_inclusion_5.45,file=paste("Res/",file_name,"_5.45_inclusion",njob,".csv",sep=""),row.names = FALSE)
write.csv(res_inclusion_opt,file=paste("Res/",file_name,"_opt_inclusion",njob,".csv",sep=""),row.names = FALSE)
write.csv(res_inclusion_opt_boot,file=paste("Res/",file_name,"_opt_boot_inclusion",njob,".csv",sep=""),row.names = FALSE)

