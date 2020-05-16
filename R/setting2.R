## simulation settings ##
## modified 2020 Mar 11 ##

df.backup<-readRDS("Dataset/bmi.cad.rds")
p<-1119 # check 
s<-1119  # check
over.dispersion<-TRUE # check

pi<-df.backup$beta.exposure
strong_ind<-order((df.backup$beta.exposure/df.backup$se.exposure)^2,decreasing = TRUE)[1:s]
strong_pi<-pi[strong_ind]
weak_ind<-1:s
weak_pi<-pi[weak_ind]
set.seed(0) # this seed is for null ses
null_se.exposure<-sample(df.backup$se.exposure,replace = TRUE,size = p-s)
sum(null_se.exposure)
null_se.outcome<-sample(df.backup$se.outcome,replace = TRUE,size=p-s)
sum(null_se.outcome)
null_se.selection<-sample(df.backup$se.selection,replace = TRUE,size=p-s)
sum(null_se.selection)
## setting 1, strong IVs ##
# data_gen<-function(){
#   useful_df<-data.frame(se.exposure=df.backup$se.exposure[strong_ind],
#                         se.outcome=df.backup$se.outcome[strong_ind],
#                         se.selection=df.backup$se.selection[strong_ind])
#   useful_df$beta.exposure<-rnorm(s,mean=strong_pi,sd = useful_df$se.exposure)
#   useful_df$beta.outcome<-rnorm(s,mean=beta0*strong_pi,sd=useful_df$se.outcome)
#   useful_df$beta.selection<-rnorm(s,mean=strong_pi,sd = useful_df$se.selection)
#   useful_df$relevant.ind<-1
#   null_df<-data.frame(se.exposure=null_se.exposure,se.outcome=null_se.outcome,se.selection=null_se.selection)
#   null_df$beta.exposure<-rnorm(p-s,0,sd=null_df$se.exposure)
#   null_df$beta.outcome<-rnorm(p-s,0,sd=null_df$se.outcome)
#   null_df$beta.selection<-rnorm(p-s,0,sd=null_df$se.selection)
#   null_df$relevant.ind<-0
#   full_df<-rbind(useful_df,null_df)
#   full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
#   full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
#   full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
#   kappa<-mean((strong_pi/useful_df$se.exposure)^2)
#   print(kappa)
#   kappa<-mean((full_df$z.exposure)^2)-1
#   print(kappa)
#   print(sum((strong_pi/useful_df$se.exposure)^2)/p)
#   #q<-Phi.tilde(lambda_our-abs(full_df$z.exposure))[1:s]+Phi.tilde(lambda_our+abs(full_df$z.exposure))[1:s]
#   #beta0*sum((full_df[1:s,])$r^2*q)
#   #sum((full_df[1:s,]$z.exposure^2*full_df[1:s,]$r^2*q))
#   return(full_df)
# }


## setting 2, many weak IV ##
# data_gen<-function(){
#   useful_df<-data.frame(se.exposure=df.backup$se.exposure[weak_ind],
#               se.outcome=df.backup$se.outcome[weak_ind],
#               se.selection=df.backup$se.selection[weak_ind])
#   print(dim(df.backup))
#   print(sum(weak_pi^2/(df.backup$se.exposure[weak_ind])^2)/p)
#   alpha_seq<-numeric(p)
#   useful_df$beta.exposure<-rnorm(s,mean=weak_pi,sd = useful_df$se.exposure)
#   useful_df$beta.outcome<-rnorm(s,mean=beta0*weak_pi+alpha_seq[1:s],sd=useful_df$se.outcome)
#   useful_df$beta.selection<-rnorm(s,mean=weak_pi,sd = useful_df$se.selection)
#   useful_df$relevant.ind<-1
#   null_df<-data.frame(se.exposure=null_se.exposure,se.outcome=null_se.outcome,se.selection=null_se.selection)
#   null_df$beta.exposure<-rnorm(p-s,0,sd=null_df$se.exposure)
#   null_df$beta.outcome<-rnorm(p-s,alpha_seq[(s+1):p],sd=null_df$se.outcome)
#   null_df$beta.selection<-rnorm(p-s,0,sd=null_df$se.selection)
#   null_df$relevant.ind<-0
#   full_df<-rbind(useful_df,null_df)
#   full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
#   full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
#   full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
#   kappa<-mean((full_df$z.exposure[1:s])^2)
#   print(kappa)
#   #q<-Phi.tilde(lambda_our-abs(full_df$z.exposure))[1:s]+Phi.tilde(lambda_our+abs(full_df$z.exposure))[1:s]
#   #beta0*sum((full_df[1:s,])$r^2*q)
#   #sum((full_df[1:s,]$z.exposure^2*full_df[1:s,]$r^2*q))
#   return(full_df)
# }

# 
# ## setting 3, strong and weak, no null ##
# data_gen<-function(){
#   full_df<-data.frame(se.exposure=df.backup$se.exposure,
#                  se.outcome=df.backup$se.outcome,
#                  se.selection=df.backup$se.selection)
#   p<-dim(df.backup)[1]
#   full_df$beta.exposure<-rnorm(p,mean=pi,sd=full_df$se.exposure)
#   full_df$beta.outcome<-rnorm(p,mean=pi*beta0,sd=full_df$se.outcome)
#   full_df$beta.selection<-rnorm(p,mean=pi,sd=full_df$se.selection)
#   full_df$relevant.ind<-1
#   full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
#   full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
#   full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
#   kappa<-mean((full_df$z.exposure)^2)-1
#   print(kappa)
#   print(sum((pi/full_df$se.exposure)^2)/p)
#   return(full_df)
# }

#all useful, with pleiotropy
data_gen<-function(){
  full_df<-data.frame(se.exposure=df.backup$se.exposure,
                      se.outcome=df.backup$se.outcome,
                      se.selection=df.backup$se.selection)
  p<-dim(df.backup)[1]
  full_df$beta.exposure<-rnorm(p,mean=pi,sd=full_df$se.exposure)
  alpha<-rnorm(p,0,2/p*sum(full_df$se.outcome))
  #alpha<-2/p*sum(full_df$se.outcome)*rlaplace(p,m=1)
  #alpha<-abs(pi)/mean(abs(pi))*rnorm(p,0,2/p*sum(full_df$se.outcome))
  full_df$beta.outcome<-rnorm(p,mean=pi*beta0+alpha,sd=full_df$se.outcome)
  full_df$beta.selection<-rnorm(p,mean=pi,sd=full_df$se.selection)
  full_df$relevant.ind<-1
  full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
  full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
  full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
  kappa<-mean((full_df$z.exposure)^2)-1
  print(kappa)
  print(sum((pi/full_df$se.exposure)^2)/p)
  return(full_df)
}

# setting 5, 200K SNPs.
# data_gen<-function(){
#   useful_df<-data.frame(se.exposure=df.backup$se.exposure[weak_ind],
#               se.outcome=df.backup$se.outcome[weak_ind],
#               se.selection=df.backup$se.selection[weak_ind])
#   print(sum(weak_pi^2/(df.backup$se.exposure[weak_ind])^2)/p)
#   alpha_seq<-numeric(p)
#   if(over.dispersion==TRUE){
#     alpha_seq<-rnorm(p,0,2/p*(sum(full_df$se.outcome)+sum(null_se.outcome)))
#   }
#   useful_df$beta.exposure<-rnorm(s,mean=weak_pi,sd = useful_df$se.exposure)
#   useful_df$beta.outcome<-rnorm(s,mean=beta0*weak_pi+alpha_seq[1:s],sd=useful_df$se.outcome)
#   useful_df$beta.selection<-rnorm(s,mean=weak_pi,sd = useful_df$se.selection)
#   useful_df$relevant.ind<-1
#   null_df<-data.frame(se.exposure=null_se.exposure,se.outcome=null_se.outcome,se.selection=null_se.selection)
#   null_df$beta.exposure<-rnorm(p-s,0,sd=null_df$se.exposure)
#   null_df$beta.outcome<-rnorm(p-s,alpha_seq[(s+1):p],sd=null_df$se.outcome)
#   null_df$beta.selection<-rnorm(p-s,0,sd=null_df$se.selection)
#   null_df$relevant.ind<-0
#   full_df<-rbind(useful_df,null_df)
#   full_df$pval.exposure<-2*pnorm(abs(full_df$beta.exposure)/full_df$se.exposure,lower.tail = FALSE)
#   full_df$pval.selection<-2*pnorm(abs(full_df$beta.selection)/full_df$se.selection,lower.tail = FALSE)
#   full_df$z.exposure<-full_df$beta.exposure/full_df$se.exposure
#   kappa<-mean((full_df$z.exposure[1:s])^2)
#   print(kappa)
#   #q<-Phi.tilde(lambda_our-abs(full_df$z.exposure))[1:s]+Phi.tilde(lambda_our+abs(full_df$z.exposure))[1:s]
#   #beta0*sum((full_df[1:s,])$r^2*q)
#   #sum((full_df[1:s,]$z.exposure^2*full_df[1:s,]$r^2*q))
#   return(full_df)
# }


