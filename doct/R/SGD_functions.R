#Gradient descent

beta_d_grad<-function(D,beta_d,Y,Npat,id,data_all_fixed,X0_inds){
  out=0 
  for (i in 1:Npat){
    inds<-which(id==i)
    cov_d<-cbind(1,Y[inds],X0_inds[inds,])
    D_ij<-D[inds]
    out=out+t(cov_d)%*%(D_ij-cov_d%*%beta_d)*(data_all_fixed$reward_i[i] - mean(data_all_fixed$reward_i))

  } 
  return( out )
}


sigma2_d_grad<-function(D,beta_d,Y,sigma2_d,Npat,id,inds_shrink_init,data_all_fixed,X0_inds){
  out=0
  for (i in 1:Npat){
    inds<-which(id==i)
    cov_d<-cbind(1,Y[inds],X0_inds[inds,])
    sum_d=(D[inds]-beta_d%*%t(cov_d) )%*%t(D[inds]-beta_d%*%t(cov_d) )
    out=out+( (  -(length(inds))/(2*sigma2_d) + sum_d/(sigma2_d^2*2) )*(data_all_fixed$reward_i[i] - mean(data_all_fixed$reward_i)) )

  }
  return(out)
}

sigma_d_grad<-function(D,beta_d,Y,sigma_d,Npat,id,inds_shrink_init,data_all_fixed,X0_inds){
  out=0
  # grad_i=numeric(Npat)
  for (i in 1:Npat){
    inds<-which(id==i)
    cov_d<-cbind(1,Y[inds],X0_inds[inds,])
    sum_d=(D[inds]-beta_d%*%t(cov_d) )%*%t(D[inds]-beta_d%*%t(cov_d) )
    out=out+( (  -(length(inds))/(sigma_d) + sum_d/(sigma_d^3) )*(data_all_fixed$reward_i[i] - mean(data_all_fixed$reward_i)) )
    
  }
  return(out)
}


data_int_func<-function(beta_v,theta_a,mu,Y,Ts,Npat,k){
  nu=exp(as.vector(t(beta_v[,1])%*%c(1,0)))
  kappa=exp(as.vector(t(beta_v[,2])%*%c(1,0)))+1
  theta_alpha=as.numeric(theta_a%*%t(cbind(1,Y)))
  alphas=k/(1+exp(theta_alpha))
  ints<-exp(mu)+alphas*dgamma(abs(diff(Ts)),shape=kappa,rate=((kappa-1)/nu) )
  
  return(ints)
}

mu_grad<-function(mu,Ts,Y,beta_v,theta_as,id,Npat,inds_shrink_init_t,end_inds,ints,data_all_fixed,ks,id_t){
  out=0
  for (i in 1:Npat){
    inds<-which(id==i)
    inds_t<-which(id_t==i)
    theta_a<-theta_as[,i]
    k<-ks[i]
    ints_new=data_int_func(beta_v,theta_a,mu,Y[inds],Ts[inds_t],Npat,k)
    out=out+( ( -exp(mu)*(Ts[end_inds[i]+i])+exp(mu)*sum(1/ints_new,na.rm=T) )*(data_all_fixed$reward_i[i] - mean(data_all_fixed$reward_i))  )
  }
  return(out) 
  
}


#' Function for calculating optimal dosage and visit parameters for a specific patient using stochastic gradient descent (SGD)
#' @param mcmc output from mcmc_joint or mcmc_separate
#' @param data_all A list of complete data. X0_inds is the baseline data for each patient, D is the dosage, Y is the longitudinal process, and Ts are the visit times,
#' id are the patient id's (ranging from 1 to the total number of patients). censor and death_time are the censoring and survival times. Npat is the total number of patients.
#' @param id_i ID of the patient to optimize
#' @param mcmc_inds Vector of post-thinning MCMC iterations indices
#' @param seed Seed.
#' @param Niter_SGD Number of SGD iterations. 
#'
#' @return A list containing the dosage and visit parameters across all SGD iterations. beta_d is the dosage linear coefficient
#' parameter vector, sigma2_d is the dosage error parameter, mu is the baseline visitation intensity, 
#' nu_1 is the visitation intensity peak, and nu_2 is the visitation intensity shape. The optimal value for beta_d is opt_beta_d 
#' and similarly for the other four parameters. 
#' 
#' @examples
#' \dontrun{
#' data("simulated_data")
#' ###########
#' #Run MCMC 
#' mcmc_settings=NULL 
#' mcmc_settings$Niter=100
#' mcmc_settings$burn.in=10
#' mcmc_settings$ndisplay=10
#' mcmc_settings$peak_dist='gamma'
#' thin=10
#' post_thin_iters<-seq(mcmc_settings$burn.in+thin,mcmc_settings$Niter,thin)
#' seed=203
#' set.seed(seed)
#' mcmc_out_joint<-mcmc_joint(simulated_data,mcmc_settings,seed)
#' ###########
#' #Run SGD
#' id_i=5
#' SGD_out<-SGD_run(mcmc_out_joint,simulated_data,id_i,post_thin_iters,seed,Niter_SGD=10)
#' 
#' }
#' @export
#' 


SGD_run<-function(mcmc,data_all,id_i,mcmc_inds,seed,Niter_SGD=1000){
  b_ils<-mcmc$b_il[,id_i,mcmc_inds]
  inds=which(data_all$id==id_i)
  init_ind<-inds[1]
  y_init=data_all$Y[init_ind]
  X0_l_i=data_all$X0_inds[init_ind,]
  dosage_bounds<-range(data_all$D)
  
  Npat<-length(mcmc_inds)
  beta_ls<-mcmc$beta_l[,mcmc_inds]
  Bs<-mcmc$B[,,mcmc_inds]
  sigma2_ls<-mcmc$sigma2_l[mcmc_inds]
  beta_ss=mcmc$beta_s[mcmc_inds]
  beta_sds=mcmc$beta_sd[mcmc_inds]
  beta_sd_cums=mcmc$beta_sd_cum[mcmc_inds]
  lambda_toxs=mcmc$lambda_tox[mcmc_inds]
  h0s=mcmc$h0[mcmc_inds]
  beta_alphas= mcmc$beta_alpha[mcmc_inds]
  shapes<-mcmc$shape[mcmc_inds]
  theta_as<-mcmc$theta_a[,mcmc_inds]
  sigma2_ds<-mcmc$sigma2_d[mcmc_inds]
  ks<-mcmc$k[mcmc_inds]
  Niter_grad<-Niter_SGD+49
  
  grad_track<-NULL
  Ncov_l<-length(X0_l_i)+4
  Ncov_d<-length(X0_l_i)+2
  Ncov_il<-3
  grad_track<-NULL 
  grad_track$id_i=id_i
  grad_track$beta_d<-array(NA, c(Ncov_d, Niter_grad))
  grad_track$sigma2_d<-rep(NA, Niter_grad)
  grad_track$tau<-rep(NA, Niter_grad)
  grad_track$beta_v<-array(0, c(2,2, Niter_grad))
  grad_track$mu<-rep(NA, Niter_grad)
  grad_track$log_mean_death<-rep(NA, Niter_grad)
  grad_track$mean_death_real<-rep(NA, Niter_grad)
  grad_track$median_death<-rep(NA, Niter_grad)
  grad_track$reward<-rep(NA, Niter_grad)
  grad_track$median_death_real<-rep(NA, Niter_grad)
  grad_track$n_visits<-rep(NA, Niter_grad)
  grad_track$learn_rate<-rep(NA, Niter_grad)
  
  grad_track$beta_d[,1]<-rowMeans(mcmc$beta_d[,mcmc_inds])
  grad_track$sigma2_d[1]<-mean(mcmc$sigma2_d[mcmc_inds])
  grad_track$tau[1]<-sqrt(grad_track$sigma2_d[1])
  # beta_d_sel<-c(1,2)
  beta_d_sel<-c(1:Ncov_d)
  beta_v<-apply(mcmc$beta_v[,,mcmc_inds],c(1,2),mean)
  mu<- mean(mcmc$mu[mcmc_inds])
  
  grad_track$beta_v[,,1]<-beta_v
  grad_track$mu[1]<-mu
  miter=2
  gradient_val<-grad_track
  
  grad_memory=50
  reward=0 
  ct=0
  beta_reward=0
  set.seed(seed)
  for (iter in miter:Niter_grad){
    if ( (iter%%50==0) & (iter>50)){
      cat("Iteration", (iter-50),'of',Niter_SGD,'complete \n')
    } 
    sigma2_ds<-rep(grad_track$sigma2_d[(iter-1)],length(mcmc_inds))
    data_all_fixed<-Data_Simu_Rcpp_single_pat_age_dgf( beta_ls,  Bs,  sigma2_ls, beta_ss,  beta_sds,
                                                       beta_sd_cums, lambda_toxs, h0s,beta_alphas,
                                                       shapes, Npat, grad_track$beta_d[,(iter-1)],
                                                       sigma2_ds, grad_track$mu[(iter-1)],
                                                       grad_track$beta_v[1,1,(iter-1)],  grad_track$beta_v[1,2,(iter-1)],
                                                       theta_as,beta_reward,b_ils, ks , X0_l_i,y_init,Ncov_l,Ncov_d,dosage_bounds)
    
    if( (!all(data_all_fixed$D>dosage_bounds[1] & data_all_fixed$D<dosage_bounds[2])) & (iter>2) ){
      print(paste0('resample iter dosage ', iter))
      grad_track$beta_d[,(iter-1)]<-grad_track$beta_d[,(iter-2)]
      data_all_fixed<-Data_Simu_Rcpp_single_pat_age_dgf( beta_ls,  Bs,  sigma2_ls, beta_ss,  beta_sds,
                                                         beta_sd_cums, lambda_toxs, h0s,beta_alphas,
                                                         shapes, Npat, grad_track$beta_d[,(iter-1)],
                                                         sigma2_ds, grad_track$mu[(iter-1)],
                                                         grad_track$beta_v[1,1,(iter-1)],  grad_track$beta_v[1,2,(iter-1)],
                                                         theta_as,beta_reward,b_ils, ks , X0_l_i,y_init,Ncov_l,Ncov_d,dosage_bounds )
    }

    D<-data_all_fixed$D
    Y<-data_all_fixed$Y
    Ts<-data_all_fixed$Ts
    id<-data_all_fixed$id
    id_t<-data_all_fixed$id_t
    X0_inds<-matrix(X0_l_i,nrow=(length(Y)),ncol=length(X0_l_i),byrow=T)
    J<-obs_count(Npat,id)
    J_cumsum<-cumsum(J) 
    J_cumsum_t<-cumsum(J+1) 
    Npat<-data_all_fixed$Npat
    EY<-data_all_fixed$EY
    ints<-1
    inds_shrink_init<-shrink_data_inds(Npat,id,init=T)
    inds_shrink_end<-shrink_data_inds(Npat,id,init=F)
    inds_shrink_init_t<-shrink_data_inds_t(Npat,id_t,init=T)
    inds_shrink_end_t<-shrink_data_inds_t(Npat,id_t,init=F)
    end_inds<-end_data_inds(Npat,id)
    init_inds<-init_data_inds(Npat,id)
    grad_track$mean_death[(iter-1)]<-mean(data_all_fixed$death_time_i)
    grad_track$log_mean_death[(iter-1)]<-mean(log(data_all_fixed$death_time_i))
    
    grad_track$median_death[(iter-1)]<-median(data_all_fixed$death_time_i)
    reward=mean(data_all_fixed$reward_i)
    grad_track$reward[iter-1]<-reward
    n_visits=length(data_all_fixed$X0)-data_all_fixed$Npat
    grad_track$n_visits[(iter-1)]=n_visits
    penalty=0.00*(n_visits)
    if (iter==2){
      init_mean=reward
      baseline=init_mean*0
      
    }
    learn_rate=10^(-2)
    init_gi=10^(-50)
    grad_track$learn_rate[(iter-1)]=learn_rate
    beta_d<-grad_track$beta_d[,(iter-1)] 
    gradient_val$beta_d[,iter]<-beta_d_grad(D,grad_track$beta_d[,(iter-1)],Y,Npat,id,data_all_fixed,X0_inds)
    beta_d_g_ii<-1/sqrt(diag(max(1,(grad_memory+1)/(iter-1))*(gradient_val$beta_d[,max(2,(iter-(grad_memory+1))):(iter-1)])%*%t(gradient_val$beta_d[,max(2,(iter-(grad_memory+1))):(iter-1)])))
    beta_d_g_ii[which(gradient_val$beta_d[,iter]==0)]=0
    if (iter<=grad_memory){
      beta_d_g_ii<-rep(init_gi,Ncov_d)
    }
    grad_track$beta_d[,iter]<-grad_track$beta_d[,(iter-1)]
    grad_track$beta_d[beta_d_sel,iter]<- grad_track$beta_d[beta_d_sel,(iter-1)]+learn_rate*(beta_d_g_ii[beta_d_sel])*gradient_val$beta_d[beta_d_sel,iter]
    
    gradient_val$tau[(iter)]<-sigma_d_grad(D,grad_track$beta_d[,(iter-1)],Y,grad_track$tau[(iter-1)],Npat,id,inds_shrink_init,data_all_fixed,X0_inds)
    # gradient_val$tau[(iter)]<-revsigma_d_grad(D,grad_track$beta_d[,(iter-1)],Y,grad_track$tau[(iter-1)],Npat,id,inds_shrink_init,data_all_fixed,X0_inds)
    tau_g_ii<-1/sqrt(diag(max(1,(grad_memory+1)/(iter-1))*t(gradient_val$tau[max(2,(iter-(grad_memory+1))):(iter-1)])%*%(gradient_val$tau[max(2,(iter-(grad_memory+1))):(iter-1)])))
    if (iter<=grad_memory){
      tau_g_ii<-init_gi
    }
    grad_track$tau[(iter)]<- (grad_track$tau[(iter-1)]+learn_rate*tau_g_ii*gradient_val$tau[(iter)])
    grad_track$sigma2_d[(iter)]<- grad_track$tau[(iter)]^2

    gradient_val$mu[(iter)]<-mu_grad(grad_track$mu[(iter-1)],Ts,Y,grad_track$beta_v[,,(iter-1)],theta_as,id,Npat,inds_shrink_init_t,end_inds,ints,data_all_fixed,ks,id_t)
    mu_g_ii<-1/sqrt(diag(max(1,(grad_memory+1)/(iter-1))*t(gradient_val$mu[max(2,(iter-(grad_memory+1))):(iter-1)])%*%(gradient_val$mu[max(2,(iter-(grad_memory+1))):(iter-1)])))
    if (iter<=grad_memory){
      mu_g_ii<-init_gi
    }
    grad_track$mu[(iter)]<-grad_track$mu[(iter-1)]+learn_rate*mu_g_ii*gradient_val$mu[(iter)]
    gradient_val$beta_v[1,2,(iter)]<-beta_v2_grad_Rcpp( J_cumsum, J_cumsum_t,  J,  Npat,  Y,  Ts, grad_track$beta_v[1,1,(iter-1)], 
                                                        grad_track$beta_v[1,2,(iter-1)],  grad_track$mu[(iter-1)],  theta_as,
                                                        ks,  data_all_fixed$reward_i) 
    
    beta_v_2_g_ii<-1/sqrt(diag(max(1,(grad_memory+1)/(iter-1))*t(gradient_val$beta_v[1,2,max(2,(iter-(grad_memory+1))):(iter-1)])%*%(gradient_val$beta_v[1,2,max(2,(iter-(grad_memory+1))):(iter-1)])))
    if (iter<=grad_memory){
      beta_v_2_g_ii<-init_gi
    } 
    
    grad_track$beta_v[1,2,(iter)]<-grad_track$beta_v[1,2,(iter-1)]+learn_rate*beta_v_2_g_ii*gradient_val$beta_v[1,2,(iter)]
    gradient_val$beta_v[1,1,(iter)]<-beta_v1_grad_Rcpp( J_cumsum, J_cumsum_t,  J,  Npat,  Y,  Ts, grad_track$beta_v[1,1,(iter-1)], 
                                                        grad_track$beta_v[1,2,(iter-1)],  grad_track$mu[(iter-1)],  theta_as,
                                                        ks,  data_all_fixed$reward_i) 
    beta_v_1_g_ii<-1/sqrt(diag(max(1,(grad_memory+1)/(iter-1))*t(gradient_val$beta_v[1,1,max(2,(iter-(grad_memory+1))):(iter-1)])%*%(gradient_val$beta_v[1,1,max(2,(iter-(grad_memory+1))):(iter-1)])))
    if (iter<=grad_memory){
      beta_v_1_g_ii<-init_gi
    }
    grad_track$beta_v[1,1,(iter)]<-grad_track$beta_v[1,1,(iter-1)]+learn_rate*beta_v_1_g_ii*gradient_val$beta_v[1,1,(iter)]
    miter<-iter
  }
  id_best<-which(grad_track$mean_reward==max(grad_track$mean_reward,na.rm=T))
  opt_beta_d<-grad_track$beta_d[,id_best]
  opt_sigma2d<-grad_track$sigma2_d[id_best]
  opt_mu<-grad_track$mu[id_best]
  opt_beta_v1<-grad_track$beta_v[1,1,id_best]
  opt_beta_v2<-grad_track$beta_v[1,2,id_best]
  
  
  output<-NULL
  output$id_i=id_i
  output$beta_d<-grad_track$beta_d[,(50:Niter_grad)]
  output$sigma2_d<-grad_track$sigma2_d[(50:Niter_grad)]
  output$mu<-grad_track$mu[(50:Niter_grad)]
  output$beta_v<-grad_track$beta_v[,,(50:Niter_grad)]
  output$nu_1<-output$beta_v[1,1,]
  output$nu_2<-output$beta_v[1,2,]
  output$opt_beta_d<-opt_beta_d
  output$opt_sigma2d<-opt_sigma2d
  output$opt_mu<-opt_mu
  output$opt_beta_v1<-opt_beta_v1
  output$opt_beta_v2<-opt_beta_v2
  
  output$median_surv<-grad_track$mean_death[(50:Niter_grad)]
  output$mean_reward<-grad_track$reward[(50:Niter_grad)]
  return(output)
}



