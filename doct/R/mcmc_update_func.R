#MCMC

sigmoid_k<-function(theta_alpha,k){
  # return(2/(1+exp(-theta_alpha)))
  return(k/(1+exp(theta_alpha)))
  
}


shrink_data_inds<-function(Npat,id,init=T){
  ind_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    if (init==T){
      l_inds<-inds[-1]
    }else{
      l_inds<-inds[-length(inds)]
    }
    ind_less<-c(ind_less,l_inds)
  }
  return(ind_less)
}

shrink_data_inds_both<-function(Npat,id){
  ind_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    l_inds<-inds[-length(inds)]
    l_inds<-l_inds[-1]

    ind_less<-c(ind_less,l_inds)
  }
  return(ind_less)
}

shrink_data_inds_t<-function(Npat,id_t,init=T){
  ind_less<-c()
  for (i in 1:Npat){
    inds<-which(id_t==i)
    if (init==T){
      l_inds<-inds[-1]
    }else{
      l_inds<-inds[-length(inds)]
    }
    ind_less<-c(ind_less,l_inds)
  }
  return(ind_less)
}

shrink_data<-function(arg,Npat,id,init=T){
  arg_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    if (init==T){
      inds_less<-inds[-1]
    }else{
      inds_less<-inds[-length(inds)]
    }
    arg_less<-c(arg_less,arg[inds_less])
  }
  return(arg_less)
}

end_data<-function(arg,Npat,id){
  arg_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    inds_less<-inds[length(inds)]
    arg_less<-c(arg_less,arg[inds_less])
  }
  return(arg_less)
}

end_data_inds<-function(Npat,id){
  arg_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    inds_less<-inds[length(inds)]
    arg_less<-c(arg_less,inds_less)
  }
  return(arg_less)
}

init_data_inds<-function(Npat,id){
  arg_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    inds_less<-inds[1]
    arg_less<-c(arg_less,inds_less)
  }
  return(arg_less)
}
inds_i<-function(Npat,id){
  sum_list<-list()
  for (i in 1:Npat){
    inds<-which(id==i)
    sum_list[[i]]<-inds
  }
  return(sum_list)
}

obs_count<-function(Npat,id){
  J<-numeric(Npat)
  for (i in 1:Npat){
    inds<-which(id==i)
    J[i]<-length(inds)
  }
  return(J)
}


data_init_death<-function(arg,death_time,Npat,id){
  arg_less<-c()
  for (i in 1:Npat){
    inds<-which(id==i)
    inds_less<-inds[-c(1,2)]
    arg_less<-c(arg_less,arg[inds_less],death_time[inds[1]])
  }
  return(arg_less)
}
loglike.surv<-function(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,
                       censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,
                       beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l){
  # alphas<-as.vector(exp(theta_a%*%t(cbind(1,Y))))
  alphas<-as.vector(sigmoid_k(theta_a%*%t(cbind(1,Y)),k)) 
  loghaz_all<-sum(log(haz_func_mixed_recent(death_time[end_inds],Ts[end_inds],beta_s,h0,EY[end_inds],shape,beta_alpha,alphas[end_inds],beta_sd,beta_sd_cum,D[end_inds],Toxicity[end_inds],lambda_tox,beta_l,b_il3[end_inds]))*censor[end_inds],na.rm=T)
  log_surv_all=loglike_surv_Rcpp( shape,beta_s,beta_alpha,h0,
                    beta_sd,beta_sd_cum, D,lambda_tox,
                    inds_shrink_end,  end_inds, EY,
                     alphas, Ts,Toxicity, death_time,Npat,length(inds_shrink_end),beta_l,b_il3, Ncov_l)
  # print(log_surv_all)
  loglike_surv<-loghaz_all+log_surv_all
  out<-NULL
  out$loglike_surv=loglike_surv
  out$log_surv_all=log_surv_all 
  return(loglike_surv)
}


loglike.vis<-function(beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist='gamma'){
  loglike_val_1<-loglike_val_2<-0
  int_all<-c()
  nu=exp(as.vector(t(beta_v[,1])%*%c(1,0)))
  kappa=exp(as.vector(t(beta_v[,2])%*%c(1,0)))+1
  a<-nu/(((kappa-1)/(kappa+1))^(1/kappa))
  alphas=sigmoid_k(as.numeric(theta_a%*%t(cbind(1,c(0,Y[-length(Y)])))),k)
  if (peak_dist=='gamma'){
  cdf_val<-pgamma(c(0,abs(diff(Ts))),shape=kappa,rate=((kappa-1)/nu))
  ints<-exp(mu)+alphas*dgamma(c(0,abs(diff(Ts))),shape=kappa,rate=((kappa-1)/nu) )
  } 
  # else if (peak_dist=='loglogistic'){
  #   cdf_val<-cdf_loglog(a,kappa,c(0,abs(diff(Ts))))
  #   ints<-exp(mu)+alphas*pdf_loglog(a,kappa,c(0,abs(diff(Ts))))
  # }

  loglike_val<-sum(log(ints[inds_shrink_init]))-sum(exp(mu)*(Ts[inds_shrink_init]-Ts[inds_shrink_end])+(alphas[inds_shrink_init])*cdf_val[inds_shrink_init])
  
  loglike_vis<-NULL
  loglike_vis$value<-loglike_val
  loglike_vis$int_all<-ints
  return(loglike_vis)
}

loglike.long<-function(Y,EY,sigma2_l,inds_shrink_init){
  loglike<-sum(dnorm(Y[inds_shrink_init],mean=EY[inds_shrink_init],sd=sqrt(sigma2_l),log=TRUE))
  return(loglike)
}

loglike.d<-function(D,beta_d,Y,X0_inds,sigma2_d){
  loglike=0 
  cov_all<-c()
  d_all<-c()
  Ncov<-length(beta_d)
  cov_all<-cbind(1,Y,X0_inds)
  loglike=sum(dnorm(D,mean=beta_d%*%t(cov_all),sd=sigma2_d,log=TRUE))
  
  return(loglike) 
}

loglike.all<-function(D,Y,EY,Ts,id,Npat,death_time,death_init_death,X0_inds,end_inds,inds_shrink_init, inds_shrink_both,
                      inds_shrink_end, sigma2_l,beta_d,sigma2_d,beta_v,theta_a,mu,k,beta_s,
                      beta_l,beta_alpha,beta_sd,beta_sd_cum,shape,h0,Toxicity,lambda_tox,b_il3,Ncov_l,censor,mean_rand){
  
  loglike_all=loglike.vis(beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k)$value+
    loglike.d(D,beta_d,Y,X0_inds,sigma2_d)+
    loglike.long(Y,EY,sigma2_l,inds_shrink_init)+
    loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,
                 censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,
                 beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  return(loglike_all)
  
}

loglike.all_separate<-function(D,Y,EY,Ts,id,Npat,death_time,death_init_death,X0_inds,end_inds,inds_shrink_init, inds_shrink_both,
                      inds_shrink_end, sigma2_l,beta_d,sigma2_d,beta_v,theta_a,mu,k,beta_s,
                      beta_l,beta_alpha,beta_sd,beta_sd_cum,shape,h0,Toxicity,lambda_tox,b_il3,Ncov_l,censor,mean_rand){
  
  loglike_all=loglike.vis(beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k)$value+
    loglike.d(D,beta_d,Y,X0_inds,sigma2_d)+
    loglike.long(Y,EY,sigma2_l,inds_shrink_init)+
    loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,Y,shape,
                 censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,
                 beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  return(loglike_all)
  
}

calculate_EY<-function(beta_l,b_il,D,Y,Ts,X0_l,Npat,inds_all){
  n<-length(Y)
  n_X0<-ncol(X0_l)
  EY<-mean_fixed<-mean_rand<-numeric(n)
  for (i in 1:Npat){
    inds<-inds_all[[i]]
    max_i<-length(inds)
    inds_less<-inds[-max_i]
    X0_l_i<-X0_l[i,]
    cov_il<-cbind(1,D[c(inds_less[1],inds_less)],Ts[inds])
    cov_l<-cbind(1,D[c(inds_less[1],inds_less)],matrix(X0_l_i,nrow=(length(inds)),ncol=n_X0,byrow=T),Ts[inds],Ts[inds]^2)
    mean_fixed[inds]<-t(beta_l)%*%t(cov_l)
    mean_rand[inds]<-t(b_il[,i])%*%t(cov_il)
    # mean_rand[inds[-1]]<-b_il[1,i]
    EY[inds]<-mean_fixed[inds]+mean_rand[inds]
  }
  out<-NULL
  out$EY<-EY
  out$fixed<-mean_fixed
  out$rand<-mean_rand
  return(out)
}

update_beta_l_separate<-function(beta_l,beta_d,beta_s,sigma2_l,sigma2_d,D,X0_l,Y,Ts,id,Npat,int_all,h0,b_il,mean_rand,B,
                        shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_both,
                        inds_shrink_init,inds_shrink_end,X0_inds,beta_sd,
                        beta_sd_cum,Toxicity,lambda_tox,b_il3,inds_all,k,Ncov_l,death_time,J_cumsum, J,accept_ind=F){
  n_X0<-ncol(X0_l)
  curr_beta_l<-beta_l
  cov_all<-cbind(1,D[inds_shrink_end],X0_inds[inds_shrink_init,],Ts[inds_shrink_init],Ts[inds_shrink_init]^2)
  Sigma<-chol2inv(chol(diag(0.01^2,Ncov_l)+t(cov_all)%*%(cov_all)/sigma2_l))
  mu<- 1/sigma2_l*Sigma%*%t(cov_all)%*%(Y[inds_shrink_init]-mean_rand[inds_shrink_init])
  prop_beta_l<-mvrnorm(n=1,mu=mu,Sigma=Sigma)
  accept_fixed<-0
  prop_out<-calculate_EY(prop_beta_l,b_il,D,Y,Ts,X0_l,Npat,inds_all)
  prop_EY<-prop_out$EY
  prop_mean_fixed<-prop_out$fixed
  curr_beta_l<-prop_beta_l
  accept_fixed<-accept_fixed+1
  curr_EY<-prop_EY
  curr_mean_fixed<-prop_out$fixed
  
  Ncov_i<-nrow(b_il)
  accept_rand<-numeric(Npat)
  
  loglikes_sum<-update_b_il( J_cumsum, J,  Npat,  D,  Y,  Ts,
                             curr_mean_fixed,  censor,  Toxicity,  b_il,  death_time,B,sigma2_l,curr_EY,theta_a,
                             beta_s, shape,  h0,  beta_alpha,  beta_sd,  beta_sd_cum,  lambda_tox,  beta_l,k, Ncov_l)
  
  b_il<-loglikes_sum$b_il
  accept_rand<-loglikes_sum$accept_rand
  
  out<-NULL
  out$beta_l<-curr_beta_l
  out$b_il<-b_il
  out$accept_fixed<-accept_fixed
  out$accept_rand<-accept_rand
  
  return(out)
}

update_beta_l<-function(beta_l,beta_d,beta_s,sigma2_l,sigma2_d,D,X0_l,Y,Ts,id,Npat,int_all,h0,b_il,mean_rand,B,
                        shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_both,
                        inds_shrink_init,inds_shrink_end,X0_inds,beta_sd,
                        beta_sd_cum,Toxicity,lambda_tox,b_il3,inds_all,k,Ncov_l,death_time,J_cumsum, J,accept_ind=F){
  n_X0<-ncol(X0_l)
  curr_beta_l<-beta_l
  cov_all<-cbind(1,D[inds_shrink_end],X0_inds[inds_shrink_init,],Ts[inds_shrink_init],Ts[inds_shrink_init]^2)
  Sigma<-chol2inv(chol(diag(0.01^2,Ncov_l)+t(cov_all)%*%(cov_all)/sigma2_l))
  mu<- 1/sigma2_l*Sigma%*%t(cov_all)%*%(Y[inds_shrink_init]-mean_rand[inds_shrink_init])
  prop_beta_l<-mvrnorm(n=1,mu=mu,Sigma=Sigma)
  accept_fixed<-0
  prop_out<-calculate_EY(prop_beta_l,b_il,D,Y,Ts,X0_l,Npat,inds_all)
  curr_out<-calculate_EY(curr_beta_l,b_il,D,Y,Ts,X0_l,Npat,inds_all)
  prop_EY<-prop_out$EY
  curr_EY<-curr_out$EY 
  curr_mean_fixed<-curr_out$fixed
  prop_mean_fixed<-prop_out$fixed
  
  loglike<-loglike.surv(beta_s,death_time,prop_beta_l,D,Y,Ts,id,Npat,h0,mean_rand,prop_EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)-
    loglike.surv(beta_s,death_time,curr_beta_l,D,Y,Ts,id,Npat,h0,mean_rand,curr_EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)

  if(runif(1)<exp(loglike)| accept_ind){
    curr_beta_l<-prop_beta_l
    accept_fixed<-accept_fixed+1
  }
  if(accept_fixed==1){
    curr_EY<-prop_EY
    curr_mean_fixed<-prop_out$fixed
  } 
  Ncov_i<-nrow(b_il)
  accept_rand<-numeric(Npat)
  
  loglikes_sum<-update_b_il( J_cumsum, J,  Npat,  D,  Y,  Ts,
               curr_mean_fixed,  censor,  Toxicity,  b_il,  death_time,B,sigma2_l,curr_EY,theta_a,
               beta_s, shape,  h0,  beta_alpha,  beta_sd,  beta_sd_cum,  lambda_tox,  beta_l,k, Ncov_l)
  
  b_il<-loglikes_sum$b_il
  accept_rand<-loglikes_sum$accept_rand
  
  out<-NULL
  out$beta_l<-curr_beta_l
  out$b_il<-b_il
  out$accept_fixed<-accept_fixed
  out$accept_rand<-accept_rand
  
  return(out)
}


update_beta_d<-function(beta_d,sigma2_d,X0_inds,Y,D,Ts,id,int_all,Npat){
  cov_all<-c()
  d_all<-c()
  Ncov<-length(beta_d)
  cov_all<-cbind(1,Y,X0_inds)
  d_all<-D
  
  Sigma<-chol2inv(chol(diag(0.01^2,Ncov)+t(cov_all)%*%(cov_all)/sigma2_d))
  mu<- 1/sigma2_d*Sigma%*%t(cov_all)%*%d_all
  new_beta_d<-mvrnorm(n=1,mu=mu,Sigma=Sigma)
  return(new_beta_d)
}

init_beta_l<-function(Ncov,sigma2_l,D,X0_l,Y,Ts,id,Npat){
  cov_all<-c()
  y_all<-c()
  n_X0<-ncol(X0_l)
  for (i in 1:Npat){
    inds<-which(id==i)
    len<-length(inds)
    d_i<-D[inds]
    X0_l_i<-X0_l[i,]
    y_i<-Y[inds]
    time<-Ts[inds]
    cov_l<-cbind(1,d_i[-len],matrix(X0_l_i,nrow=(length(inds)-1),ncol=n_X0,byrow=T),time[-1],time[-1]^2)
    cov_all<-rbind(cov_all,cov_l)
    y_all<-c(y_all,y_i[-1])
  }
  Sigma<-chol2inv(chol(diag(0.01^2,Ncov)+t(cov_all)%*%(cov_all)/sigma2_l))
  mu<- 1/sigma2_l*Sigma%*%t(cov_all)%*%y_all
  return(mvrnorm(n=1,mu=mu,Sigma=Sigma))
}


update_beta_s<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_s,h0,mean_rand,EY,
                        shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,
                        beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l){
  beta_s_out<-beta_s
  beta_s_new<-rnorm(1,beta_s,jump_s)
  accept=0
  prop_loglik<-loglike.surv(beta_s_new,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,
                            beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  
  if (is.na(curr_loglik_surv)){ 
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,
                                   beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dnorm(beta_s_new,mean=0,sd=100,log=T)-dnorm(beta_s,mean=0,sd=100,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dnorm(beta_s_new,mean=0,sd=100,log=T)-dnorm(beta_s,mean=0,sd=100,log=T)+ prop_loglik-curr_loglik_surv
  }
  if (runif(1)<exp(loglike)){
    beta_s_out<-beta_s_new
    accept=accept+1
    curr_loglik_surv<-prop_loglik
  } 
  out<-NULL
  out$beta_s<-beta_s_out 
  out$accept<-accept 
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_lambda_tox<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_lambda_tox,h0,mean_rand,
                            EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,
                            beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l,sd_prior=NA){
  lambda_tox_out<-lambda_tox
  Tox_out<-Toxicity
  lambda_tox_new<-rnorm(1,lambda_tox,jump_lambda_tox)

  if (lambda_tox_new<1){
    # lambda_tox_new=abs(lambda_tox_new)
    lambda_tox_new=lambda_tox
  }
  prior_shape=0.001
  prior_rate=0.001
  
  accept=0
  Toxicity_new<-toxicity_calc_all(D,Ts,id,lambda_tox_new,Npat)
  prop_loglik<- loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,
                             theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,
                             beta_sd_cum,Toxicity_new,lambda_tox_new,b_il3,k,Ncov_l)
  if (is.na(curr_loglik_surv)){
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,
                                   beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dgamma(lambda_tox_new,shape=prior_shape,rate=prior_rate,log=T)-dgamma(lambda_tox,shape=prior_shape,rate=prior_rate,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dgamma(lambda_tox_new,shape=prior_shape,rate=prior_rate,log=T)-dgamma(lambda_tox,shape=prior_shape,rate=prior_rate,log=T)+ prop_loglik-curr_loglik_surv
  }
  
  if (runif(1)<exp(loglike)){
    lambda_tox_out<-lambda_tox_new 
    Tox_out<-Toxicity_new 
    accept=accept+1
    curr_loglik_surv<-prop_loglik
    
  } 
  out<-NULL
  out$lambda_tox<-lambda_tox_out 
  out$Toxicity<-Tox_out 
  out$accept<-accept 
  out$curr_loglik_surv<-curr_loglik_surv
  
  return(out)
}
update_beta_sd<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_sd,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l){
  beta_sd_out<-beta_sd
  beta_sd_new<-rnorm(1,beta_sd,jump_sd) 
  accept=0
  prop_loglik<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,
                            beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd_new,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  if (is.na(curr_loglik_surv)){
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,
                                   end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dnorm(beta_sd_new,mean=0,sd=100,log=T)-dnorm(beta_sd,mean=0,sd=100,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dnorm(beta_sd_new,mean=0,sd=100,log=T)-dnorm(beta_sd,mean=0,sd=100,log=T)+ prop_loglik-curr_loglik_surv
  }
  
  if (runif(1)<exp(loglike)){
    beta_sd_out<-beta_sd_new 
    accept=accept+1
    curr_loglik_surv<-prop_loglik
  }  
  out<-NULL
  out$beta_sd<-beta_sd_out
  out$accept<-accept  
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_beta_sd_cum<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_sd_cum,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l){
  beta_sd_cum_out<-beta_sd_cum
  beta_sd_cum_new<-rnorm(1,beta_sd_cum,jump_sd_cum)
  accept=0
  prop_loglik<- loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,
                             end_inds,inds_shrink_end, beta_sd,beta_sd_cum_new,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  
  if (is.na(curr_loglik_surv)){
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,
                                   end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dnorm(beta_sd_cum_new,mean=0,sd=100,log=T)-dnorm(beta_sd_cum,mean=0,sd=100,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dnorm(beta_sd_cum_new,mean=0,sd=100,log=T)-dnorm(beta_sd_cum,mean=0,sd=100,log=T)+ prop_loglik-curr_loglik_surv
  }
  if (runif(1)<exp(loglike)){
    beta_sd_cum_out<-beta_sd_cum_new
    accept=accept+1
    curr_loglik_surv<-prop_loglik
  }
  out<-NULL
  out$beta_sd_cum<-beta_sd_cum_out
  out$accept<-accept
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_h0<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_h0,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l){
  h0_out<-h0
  h0_new<-rnorm(1,h0,jump_h0)
  accept=0
  prop_loglik<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0_new,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,
                            end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
 
  if (is.na(curr_loglik_surv)){
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,
                                   end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dnorm(h0_new,mean=0,sd=100,log=T)-dnorm(h0_out,mean=0,sd=100,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dnorm(h0_new,mean=0,sd=100,log=T)-dnorm(h0_out,mean=0,sd=100,log=T)+ prop_loglik-curr_loglik_surv
  }
  if (runif(1)<exp(loglike)){
    h0_out<-h0_new
    accept=accept+1
    curr_loglik_surv<-prop_loglik
  }
  out<-NULL
  out$h0<-h0_out
  out$accept<-accept
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_shape<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_shape,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l){
  shape_out<-shape
  shape_new<-rnorm(1,shape,jump_shape)
  if(shape_new<0){
    # shape_new= -shape_new
    shape_new= shape
    
  }
  prior_shape=0.01
  prior_rate=0.01
  accept=0
  prop_loglik<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape_new,censor,theta_a,
                            beta_alpha,death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  
  if (is.na(curr_loglik_surv)){
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,
                                   death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dgamma(shape_new,shape=prior_shape,rate=prior_rate,log=T)-dgamma(shape_out,shape=prior_shape,rate=prior_rate,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dgamma(shape_new,shape=prior_shape,rate=prior_rate,log=T)-dgamma(shape_out,shape=prior_shape,rate=prior_rate,log=T)+ prop_loglik-curr_loglik_surv
  }
  
  if (runif(1)<exp(loglike)){
    shape_out<-shape_new
    accept=accept+1
    curr_loglik_surv<-prop_loglik
  }
  out<-NULL
  out$shape<-shape_out
  out$shape_new<-shape_new
  out$accept<-accept
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_beta_alpha<-function(beta_s,beta_l,D,Y,Ts,id,Npat,death_time,jump_beta_alpha,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,end_inds,inds_shrink_end,beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,curr_loglik_surv=NA,k,Ncov_l){
  beta_alpha_out<-beta_alpha
  beta_alpha_new<-rnorm(1,beta_alpha,jump_beta_alpha)
  # beta_alpha_new<-max(0,rnorm(1,beta_alpha,jump_beta_alpha))
  accept=0
  prop_loglik<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha_new,
                            death_init_death,end_inds,inds_shrink_end,beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
  
  if (is.na(curr_loglik_surv)){
    curr_loglik_surv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,
                                   death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<-dnorm(beta_alpha_new,mean=0,sd=100,log=T)-dnorm(beta_alpha_out,mean=0,sd=100,log=T)+prop_loglik-curr_loglik_surv
  } else{
    loglike<-dnorm(beta_alpha_new,mean=0,sd=100,log=T)-dnorm(beta_alpha_out,mean=0,sd=100,log=T)+ prop_loglik-curr_loglik_surv
  }
  if (runif(1)<exp(loglike)){
    beta_alpha_out<-beta_alpha_new
    accept=accept+1
    curr_loglik_surv<-prop_loglik
  }
  out<-NULL
  out$beta_alpha<-beta_alpha_out
  # out$beta_alpha_new<-beta_alpha_new
  out$accept<-accept
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_beta_v<-function(beta_v,Y,Ts,mu,theta_a,Npat,jump_v,id,inds_shrink_init,inds_shrink_end,k,peak_dist){
  cur_beta_v<-beta_v
  accept<-matrix(c(0,0,0,0),nrow=2)
  for(j in 1:2){
    for (i in 1){
      beta_v_orig<-cur_beta_v[i,j]
      prop_beta_v<-cur_beta_v 
      prop_beta_v[i,j]<-rnorm(1,beta_v_orig,jump_v[i,j])
      loglike<- dnorm(prop_beta_v[i,j],mean=0,sd=100,log=T)-dnorm(beta_v_orig,mean=0,sd=100,log=T)+
        loglike.vis(prop_beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist)$value-
        loglike.vis(cur_beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist)$value
      # loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,prop_beta_v,beta_nu)-
      # loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,cur_beta_v,beta_nu)
      
      if (runif(1)<exp(loglike)){
        cur_beta_v<-prop_beta_v
        accept[i,j]<-accept[i,j]+1
      }
    }
  }
  out<-NULL
  out$beta_v<-cur_beta_v
  out$accept<-accept
  return(out)
}

update_mu<-function(beta_v,Y,Ts,mu,theta_a,Npat,jump_mu,id,inds_shrink_init,inds_shrink_end,k,peak_dist){
  cur_mu<-mu
  accept<-0
  prop_mu<-rnorm(1,cur_mu,jump_mu)
  loglike<- dnorm(prop_mu,mean=0,sd=100,log=T)-dnorm(cur_mu,mean=0,sd=100,log=T)+
    loglike.vis(beta_v,theta_a,prop_mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist)$value-
    loglike.vis(beta_v,theta_a,cur_mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist)$value
  if (runif(1)<exp(loglike) & prop_mu<0){
    cur_mu<-prop_mu
    accept<-accept+1
  }
  out<-NULL
  out$mu<-cur_mu
  out$accept<-accept
  return(out)
} 


update_theta_a<-function(beta_v,Y,Ts,mu,theta_a,Npat,jump_theta_a,beta_alpha,beta_s,death_time,beta_l,D,id,h0,mean_rand,EY,shape,
                         censor,inds_shrink_init,inds_shrink_end,death_init_death,end_inds,beta_sd,beta_sd_cum,
                         Toxicity,lambda_tox,inds_shrink_both,b_il3,curr_loglik_surv=NA,k,Ncov_l,peak_dist){
  cur_theta_a<-theta_a
  accept<-c(0,0)
  for (i in 1:2){
    theta_a_orig<-cur_theta_a[i]
    prop_theta_a<-cur_theta_a 
    prop_theta_a[i]<-rnorm(1,theta_a_orig,jump_theta_a[i])
    prop_logsurv<-loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,prop_theta_a,beta_alpha,
                               death_init_death,end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,k,Ncov_l)
    loglike<- dnorm(prop_theta_a[i],mean=0,sd=100,log=T)-dnorm(theta_a_orig,mean=0,sd=100,log=T)+
      loglike.vis(beta_v,prop_theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist)$value-
      loglike.vis(beta_v,cur_theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,k,peak_dist)$value+
      prop_logsurv-curr_loglik_surv     
    if (runif(1)<exp(loglike)){ 
      cur_theta_a<-prop_theta_a
      accept[i]<-accept[i]+1
      curr_loglik_surv<-prop_logsurv
    }
  }
  out<-NULL
  out$theta_a<-cur_theta_a
  out$accept<-accept
  out$curr_loglik_surv<-curr_loglik_surv
  return(out)
}

update_k<-function(beta_v,Y,Ts,mu,theta_a,Npat,jump_theta_a,beta_alpha,beta_s,death_time,beta_l,D,id,h0,mean_rand,EY,shape,
                         censor,inds_shrink_init,inds_shrink_end,death_init_death,end_inds,beta_sd,beta_sd_cum,Toxicity,
                   lambda_tox,inds_shrink_both,b_il3,curr_loglik_surv=NA,k,jump_k,a=400,b=200,Ncov_l,peak_dist){
  cur_k<-k
  accept<-0
  prop_k<-rnorm(1,cur_k,jump_k)

  if (prop_k<0){
    prop_k=k
  }
  loglike<- dgamma(prop_k, shape=a,rate=b,log=T)-dgamma(cur_k, shape=a,rate=b,log=T)+
    loglike.vis(beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,prop_k,peak_dist)$value-
    loglike.vis(beta_v,theta_a,mu,Y,Ts,id,Npat,inds_shrink_init,inds_shrink_end,cur_k,peak_dist)$value+
    loglike.surv(beta_s,death_time,beta_l,D,Y,Ts,id,Npat,h0,mean_rand,EY,shape,censor,theta_a,beta_alpha,death_init_death,
                 end_inds,inds_shrink_end, beta_sd,beta_sd_cum,Toxicity,lambda_tox,b_il3,prop_k,Ncov_l)-
    curr_loglik_surv     
  if (runif(1)<exp(loglike)){ 
    cur_k<-prop_k
    accept<-accept+1
  }
  out<-NULL
  out$k<-cur_k
  out$accept<-accept 
  return(out)
}


update_sigma2_l<-function(lambda1_l,lambda2_l,mean_fixed,Y,Npat,id,B,inds_shrink_init){
  n<-length(Y)-Npat
  sum=0
  sum=(Y[inds_shrink_init]-mean_fixed[inds_shrink_init])%*%(Y[inds_shrink_init]-mean_fixed[inds_shrink_init])/2
  lambda1_l_post<-lambda1_l+n/2
  lambda2_l_post<-lambda2_l+sum
  new_sigma2_l<-1/rgamma(1,shape=lambda1_l_post,rate=lambda2_l_post)
  return(new_sigma2_l) 
}

update_B<-function(c_B,d_B,b_il,Npat){
  n<-Npat
  # C_post<-c_B+n
  # D_post<-d_B+b_il%*%t(b_il)
  C_post<-n
  D_post<-b_il%*%t(b_il)
  new_B<- rinvwishart(C_post,(D_post))

  return(new_B)
}

update_sigma2_d<-function(lambda1_d,lambda2_d,Y,Npat,id,int_all,beta_d,B,X0_d_inds,D){
  n<-length(Y)-Npat
  cov_d<-cbind(1,Y,X0_d_inds)
  
  sum=(D-beta_d%*%t(cov_d) )%*%t(D-beta_d%*%t(cov_d) )/2
  lambda1_d_post<-lambda1_d+n/2
  lambda2_d_post<-lambda2_d+sum
  new_sigma2_d<-1/rgamma(1,shape=lambda1_d_post,rate=lambda2_d_post)
  return(new_sigma2_d)
}



