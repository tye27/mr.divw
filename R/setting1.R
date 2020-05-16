## simulation settings ##
## modified 2020 Mar 11 ##
#setwd("~/Dropbox (Penn)/!! My experience/Research/IV/Simulation/2020 Mar revision")
######################################
n<-1e4
p<-2e3
s<-2e3  # check, set as 50, 1000, 4000
true_var<-TRUE # check, know the true variance
over.dispersion<-FALSE # check  
h<-sqrt(0.2) # 20% total heritability 
var_factor<-1
gamma_coef<-readRDS(file="Dataset/standard_normal_random_effect.rds")[1:s]
gamma_coef<-gamma_coef*(h/sqrt(s))*sqrt(2)
#gamma_coef[1:s]<-0.5
######################################

beta0<-0.4
if(s<p){
  gamma_coef[(s+1):p]<-0  
}


var_x<-0.5*sum(gamma_coef[1:s]^2 )+(1-h^2)*var_factor
var_y<-beta0^2*var_x+(2*beta0+1)*(1-h^2)*0.6*var_factor+var_factor
sigma_yj_true<-sqrt((n*0.5)^(-1)*(var_y-beta0^2*gamma_coef[1:p]^2*0.5))
sigma_xj_true<-sqrt((n*0.5)^(-1)*(var_x-gamma_coef[1:p]^2*0.5))
kappa<-sum( gamma_coef[1:p]^2/ sigma_xj_true^2)/p
round(kappa,3)
round(kappa*sqrt(p),3)
round(kappa*p/n,3)

file_name<-paste0("individual_", true_var, "n=",n, "p=",p, "s=",s, "h^2=", h^2, "beta0=",beta0)
data_gen_onesample<-function(){
  z<-matrix(nrow=n,ncol=p)
  for(j in 1:p){
    tmp<-rmultinom(n,1,c(0.25,0.5,0.25)) # z's are independent
    z[,j]<-0*(tmp[1,]==1)+1*(tmp[2,]==1)+2*(tmp[3,]==1)
  }
  u<-rnorm(n,0,sqrt((1-h^2)*0.6))
  x<-z[,1:s] %*% gamma_coef[1:s] +u+rnorm(n,0,sqrt((1-h^2)*0.4))
  y<-10+beta0*x+u+rnorm(n,0,sqrt(var_factor))
  return(list(z=z,x=x,y=y))
}


data_gen<-function(true_var=FALSE){
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

# data_gen_indep<-function(true_var){
#   full_df<-data.frame(beta.exposure=numeric(p))
#   full_df$beta.exposure<-rnorm(p,gamma_coef, sigma_xj_true)
#   full_df$beta.outcome<-rnorm(p,beta0*gamma_coef, sigma_yj_true)
#   if(true_var==TRUE){
#     full_df$se.exposure<-sigma_xj_true
#     full_df$se.outcome<-sigma_yj_true
#   }else if (true_var==FALSE){
#     tmp1<-data_gen_onesample() # exposure dataset
#     tmp2<-data_gen_onesample() # outcome dataset
#     full_df$se.exposure<-sapply(1:p, function(i) {
#       se <- summary( lm(tmp1$x ~ tmp1$z[,i]))$coefficients[2,2]
#       return(se)
#     })
#     full_df$se.outcome<-sapply(1:p, function(i) {
#       se <- summary( lm(tmp2$y ~ tmp2$z[,i]))$coefficients[2,2]
#       return(se)
#     })
#   }
#   return(full_df)
# }
