#' MCMC function for joint longitudinal, survival, and dosage and visitation model
#'
#' @param data_all A list of complete data. X0_inds is the baseline data for each patient, D is the dosage, Y is the longitudinal process, and Ts are the visit times,
#' id are the patient id's (ranging from 1 to the total number of patients), censor and surv_time are the censoring and survival times, and Npat is the total number of patients.
#' @param mcmc_settings A list for MCMC setup. burn.in is the number of total iterations to burn, ndisplay is number of iterations per which the display message will
#' appear, and Niter is the total number of iterations (including burn-in iterations).
#' @param seed Seed for running MCMC 
#'
#'
#' @return 
#'   \item{beta_l}{Linear coefficient parameter in longitudinal submodel}
#'   \item{sigma2_l}{Error parameter in longitudinal submodel}
#'   \item{beta_d}{Linear coefficient parameter in dosing submodel}
#'   \item{sigma2_d}{Error parameter in dosing submodel}
#'   \item{nu_1}{Visitation intensity peak parameter in vistation submodel}
#'   \item{nu_2}{Visitation intensity shape parameter in vistation submodel}
#'   \item{mu}{Baseline visitation intensity parameter in vistation submodel}
#'   \item{beta_s1}{Creatinine-associated parameter in survival submodel}
#'   \item{beta_s2}{Dosage-associated parameter in survival submodel}
#'   \item{beta_s3}{Toxicity-associated parameter in survival submodel}
#'   \item{beta_s4}{Visitation-associated parameter in survival submodel}
#'   \item{h0}{Baseline hazard parameter in survival submodel}
#'   \item{shape}{Shape parameter in survival submodel}
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
#' 
#' }
#' 
#' @export
#'  

mcmc_joint<-function(data_all,mcmc_settings,seed){
  Niter=mcmc_settings$Niter
  burn.in=mcmc_settings$burn.in
  ndisplay=mcmc_settings$ndisplay
  peak_dist=mcmc_settings$peak_dist
  D<-data_all$D
  Y<-data_all$Y
  Ts<-data_all$Ts
  id<-data_all$id
  censor<-data_all$censor
  surv_time<-data_all$surv_time
  Npat<-data_all$Npat
  
  #Prepping data
  inds_shrink_init<-shrink_data_inds(Npat,id,init=T)
  inds_shrink_end<-shrink_data_inds(Npat,id,init=F)
  inds_shrink_both<-shrink_data_inds_both(Npat,id)
  end_inds<-end_data_inds(Npat,id)
  inds_all<-inds_i(Npat,id)
  J<-obs_count(Npat,id)
  J_cumsum<-cumsum(J) 
  death_init_death<-data_init_death(Ts,surv_time,Npat,id)
  X0_l<-data_all$X0_inds[end_inds,]
  X0_d<-data_all$X0_inds[end_inds,]
  X0_inds<-data_all$X0_inds
  X0<-X<-X0_l
  Ncov_l<-ncol(X0_l)+4 
  Ncov_d<-ncol(X0_d)+2 
  Ncov_il<-3
  mcmc<-NULL 
  mcmc$Ncov_l=Ncov_l
  mcmc$Ncov_d=Ncov_d
  mcmc$Ncov_il=Ncov_il
  mcmc$beta_l<-array(NA, c(Ncov_l, Niter))
  mcmc$b_il<-array(NA, c(Ncov_il,Npat,Niter))
  mcmc$B<-array(NA,c(Ncov_il,Ncov_il,Niter))
  mcmc$beta_d<-array(NA, c(Ncov_d, Niter))
  mcmc$sigma2_l<-rep(NA, Niter)
  mcmc$sigma2_d<-rep(NA, Niter)
  mcmc$lambda_tox<-rep(NA, Niter)
  mcmc$k<-rep(NA,Niter)
  mcmc$beta_s<-rep(NA, Niter)
  mcmc$beta_sd<-rep(NA, Niter)
  mcmc$beta_sd_cum<-rep(NA, Niter)
  mcmc$beta_v<-array(NA, c(2,2, Niter))
  mcmc$theta_a<-array(NA, c(2,Niter))
  mcmc$mu<-rep(NA, Niter)
  mcmc$h0<-rep(NA, Niter)
  mcmc$shape<-rep(NA, Niter)
  mcmc$beta_alpha<-rep(NA, Niter)

  mcmc$beta_l[,1]<-init_beta_l(Ncov_l,0.5^2,D,X0_l,Y,Ts,id,Npat)
  mcmc$b_il[,,1]<-matrix(0,nrow=Ncov_il,ncol=Npat)
  b_il3<-mcmc$b_il[3,id,1]
  mcmc$lambda_tox[1]<-60
  mcmc$k[1]<-2
  mcmc$beta_s[1]<- 0
  mcmc$beta_sd[1]<- 0
  mcmc$beta_sd_cum[1]<- 0
  mcmc$h0[1]<- 0
  mcmc$shape[1]<- 1
  mcmc$beta_alpha[1]<- 0
  mcmc$beta_v[,,1]<-matrix(c(2,0,2,0),nrow=2)
  mcmc$beta_d[,1]<-0
  mcmc$beta_d[1,1]<-5
  mcmc$sigma2_l[1]<-0.5^2
  mcmc$sigma2_d[1]<-0.5^2
  mcmc$B[,,1]<-matrix(c(0.3^2,0,0,
                        0,0.2^2,0,
                        0,0,(5*10^(-4))^2),nrow=3)
  mcmc$theta_a[,1]<- c(9,-1.4)
  mcmc$mu[1]<- -4.5
  
  # mcmc$theta_a[,1]<- c(11.5,-2.05)
  # mcmc$mu[1]<- -4
  
  #hyperparameter values
  lambda1_l<-lambda2_l<-lambda1_d<-lambda2_d<-0.01
  a<-400
  b<-200
  c_B<-d_B<-0
  #MH settings + acceptance rates
  accept_beta_l<-0
  accept_beta_v<-matrix(c(0,0,0,0),nrow=2)
  accept_beta_s<-0 
  accept_beta_sd<-0
  accept_lambda_tox<-0
  accept_beta_sd_cum<-0 
  accept_mu<-0
  accept_theta_a<-c(0,0) 
  accept_h0<-0
  accept_shape<-0 
  accept_beta_alpha<-0
  accept_k<-0
  jump_v<-matrix(c(0.07,0,0.8,0),nrow=2)/3 
  jump_s<-0.06  
  jump_sd<-0.18
  jump_sd_cum<-0.18
  jump_lambda_tox<-7
  jump_mu<-0.05
  jump_theta_a<-c(0.4,0.08)/4.5
  jump_h0<-0.25 
  jump_k<-0.05 
  jump_shape<-0.04
  jump_beta_alpha<-0.4 
  miter=2  

  set.seed(seed)
  for (iter in miter:Niter){
    if (iter%%ndisplay==0){
      cat("Iteration", iter,'of',Niter,'complete \n') 
      # print(mcmc$lambda_tox[iter-1])
      
    } 
    EY_out<-calculate_EY(mcmc$beta_l[,iter-1],mcmc$b_il[,,iter-1],D,Y,Ts,X0_l,Npat,inds_all)
    EY<-EY_out$EY
    mean_rand<-EY_out$rand
    int_all<-NA
    # beta_l<-mcmc$beta_l[,iter-1] ; beta_d<- mcmc$beta_d[,iter-1] ; beta_s<-mcmc$beta_s[iter-1] ; sigma2_l<-mcmc$sigma2_l[iter-1] ;
    # sigma2_d<-mcmc$sigma2_d[iter-1] ;  h0<- mcmc$h0[iter-1] ; b_il<- mcmc$b_il[,,iter-1] ; B<-mcmc$B[,,iter-1]
    # beta_s<-mcmc$beta_s[iter-1] ; h0<-mcmc$h0[iter-1]; shape<-mcmc$shape[iter-1];
    # theta_a<-mcmc$theta_a[,(iter-1)];beta_alpha<-mcmc$beta_alpha[iter-1];beta_sd=mcmc$beta_sd[iter-1];
    # beta_sd_cum=mcmc$beta_sd_cum[iter-1];Toxicity<-toxicity_calc_all(D,Ts,id,mcmc$lambda_tox[iter-1] ) ;lambda_tox=mcmc$lambda_tox[iter-1] ; k<-mcmc$k[iter-1]
    # k<-mcmc$k[iter-1]
    Toxicity<-toxicity_calc_all(D,Ts,id,mcmc$lambda_tox[iter-1],Npat )
    out<-update_beta_l(mcmc$beta_l[,iter-1],mcmc$beta_d[,iter-1],mcmc$beta_s[iter-1],mcmc$sigma2_l[iter-1],
                       mcmc$sigma2_d[iter-1],D,X0_l,Y,Ts,id,Npat,int_all=NA,mcmc$h0[iter-1],mcmc$b_il[,,iter-1],
                       mean_rand,mcmc$B[,,iter-1],mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],
                       death_init_death,end_inds,inds_shrink_both, inds_shrink_init,inds_shrink_end,X0_inds,
                       mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],
                       Toxicity,mcmc$lambda_tox[iter-1],b_il3,inds_all,mcmc$k[iter-1],Ncov_l,surv_time,J_cumsum, J)
    mcmc$beta_l[,iter]<-out$beta_l
    mcmc$b_il[,,iter]<-out$b_il
    b_il3<-mcmc$b_il[3,id,iter]
    
    mcmc$B[,,iter]<-update_B(c_B,d_B,mcmc$b_il[,,iter],Npat)
    EY_out<-calculate_EY(mcmc$beta_l[,iter],mcmc$b_il[,,iter],D,Y,Ts,X0_l,Npat,inds_all)
    EY<-EY_out$EY
    mean_rand<-EY_out$rand
    mean_fixed<-EY_out$fixed
    # 
    # beta_s<-mcmc$beta_s[iter-1] ; beta_l<-mcmc$beta_l[,iter]; h0<-mcmc$h0[iter-1]; shape<-mcmc$shape[iter-1]; theta_a<-mcmc$theta_a[,(iter-1)];beta_alpha<-mcmc$beta_alpha[iter-1];beta_sd=mcmc$beta_sd[iter-1]; beta_sd_cum=mcmc$beta_sd_cum[iter-1];Toxicity<-toxicity_calc_all(D,Ts,id,mcmc$lambda_tox[iter-1] ) ;lambda_tox=mcmc$lambda_tox[iter-1]
    #  
    curr_loglik_surv=NA
    out_lambda_tox<-update_lambda_tox(mcmc$beta_s[iter-1],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,
                                      jump_lambda_tox,mcmc$h0[iter-1],mean_rand,EY,mcmc$shape[iter-1],
                                      censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,
                                      inds_shrink_end,mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],Toxicity,mcmc$lambda_tox[iter-1],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$lambda_tox[iter]<-out_lambda_tox$lambda_tox
    jump_lambda_tox<- sqrt(out_lambda_tox$lambda_tox)*5
    Toxicity=out_lambda_tox$Toxicity
    curr_loglik_surv<-out_lambda_tox$curr_loglik_surv
    
    out_s<-update_beta_s(mcmc$beta_s[iter-1],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_s,mcmc$h0[iter-1],mean_rand,EY,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_s[iter]<-out_s$beta_s
    curr_loglik_surv<-out_s$curr_loglik_surv
    
    out_sd<-update_beta_sd(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_sd,mcmc$h0[iter-1],mean_rand,EY,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_sd[iter]<-out_sd$beta_sd
    curr_loglik_surv<-out_sd$curr_loglik_surv
    
    out_sd_cum<-update_beta_sd_cum(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_sd_cum,mcmc$h0[iter-1],mean_rand,EY,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter-1],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_sd_cum[iter]<-out_sd_cum$beta_sd_cum
    curr_loglik_surv<-out_sd_cum$curr_loglik_surv
    
    out_h0<-update_h0(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_h0,mcmc$h0[iter-1],mean_rand,EY,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$h0[iter]<-out_h0$h0
    curr_loglik_surv<-out_h0$curr_loglik_surv
    
    out_shape<-update_shape(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_shape,mcmc$h0[iter],mean_rand,EY,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$shape[iter]<-out_shape$shape
    curr_loglik_surv<-out_shape$curr_loglik_surv
    
    out_beta_alpha<-update_beta_alpha(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_beta_alpha,mcmc$h0[iter],mean_rand,EY,mcmc$shape[iter],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_alpha[iter]<-out_beta_alpha$beta_alpha
    curr_loglik_surv<-out_beta_alpha$curr_loglik_surv
    
    out_v<-update_beta_v(mcmc$beta_v[,,iter-1],Y,Ts,mcmc$mu[iter-1],mcmc$theta_a[,iter-1],Npat,jump_v,id,inds_shrink_init,inds_shrink_end,mcmc$k[iter-1],peak_dist)
    mcmc$beta_v[,,iter]<-out_v$beta_v
    
    out_mu<-update_mu(mcmc$beta_v[,,iter],Y,Ts,mcmc$mu[iter-1],mcmc$theta_a[,iter-1],Npat,jump_mu,id,inds_shrink_init,inds_shrink_end,mcmc$k[iter-1],peak_dist)
    mcmc$mu[iter]<-out_mu$mu
    
    out_theta_a<-update_theta_a(mcmc$beta_v[,,iter],Y,Ts,mcmc$mu[iter],mcmc$theta_a[,iter-1],Npat,jump_theta_a,mcmc$beta_alpha[iter],mcmc$beta_s[iter],
                                surv_time,mcmc$beta_l[,iter],D,id,mcmc$h0[iter],mean_rand,EY,mcmc$shape[iter],censor,inds_shrink_init,inds_shrink_end,death_init_death,end_inds,
                                mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity,mcmc$lambda_tox[iter],inds_shrink_both,b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l,peak_dist)
    mcmc$theta_a[,iter]<-out_theta_a$theta_a
    curr_loglik_surv<-out_theta_a$curr_loglik_surv
    
    out_k<-update_k(mcmc$beta_v[,,iter],Y,Ts,mcmc$mu[iter],mcmc$theta_a[,iter],Npat,jump_theta_a,mcmc$beta_alpha[iter],mcmc$beta_s[iter],
                    surv_time,mcmc$beta_l[,iter],D,id,mcmc$h0[iter],mean_rand,EY,mcmc$shape[iter],censor,inds_shrink_init,inds_shrink_end,death_init_death,end_inds,
                    mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity,mcmc$lambda_tox[iter],inds_shrink_both,b_il3,curr_loglik_surv,mcmc$k[iter-1],jump_k,a,b,Ncov_l,peak_dist)
    mcmc$k[iter]<-out_k$k
    
    mcmc$beta_d[,iter]<-update_beta_d(mcmc$beta_d[,iter-1],mcmc$sigma2_d[iter-1],X0_inds,Y,D,Ts,id,int_all,Npat)
    
    mcmc$sigma2_l[iter]<-update_sigma2_l(lambda1_l,lambda2_l,EY,Y,Npat,id,mcmc$B[iter],inds_shrink_init)
    mcmc$sigma2_d[iter]<-update_sigma2_d(lambda1_d,lambda2_d,Y,Npat,id,int_all,mcmc$beta_d[,iter],mcmc$B[iter],X0_inds,D)
    
    if (iter>burn.in){
      accept_beta_l<-accept_beta_l+out$accept_fixed
      accept_beta_v<-accept_beta_v+out_v$accept 
      accept_beta_s<-accept_beta_s+out_s$accept
      accept_lambda_tox<-accept_lambda_tox+out_lambda_tox$accept
      
      accept_beta_sd<-accept_beta_sd+out_sd$accept 
      accept_beta_sd_cum<-accept_beta_sd_cum+out_sd_cum$accept
      accept_h0<-accept_h0+out_h0$accept
      accept_shape<-accept_shape+out_shape$accept
      accept_beta_alpha<-accept_beta_alpha+out_beta_alpha$accept
      accept_mu<-accept_mu+out_mu$accept
      accept_k<-accept_k+out_k$accept
      accept_theta_a<-accept_theta_a+out_theta_a$accept
    }
    
    miter<-iter
  }
  mcmc$accept_beta_l<-accept_beta_l
  mcmc$accept_beta_v<-accept_beta_v
  mcmc$accept_beta_s<-accept_beta_s
  mcmc$accept_lambda_tox<-accept_lambda_tox
  mcmc$accept_beta_sd<-accept_beta_sd
  mcmc$accept_beta_sd_cum<-accept_beta_sd_cum
  mcmc$accept_h0<-accept_h0
  mcmc$accept_shape<-accept_shape
  mcmc$accept_beta_alpha<-accept_beta_alpha
  mcmc$accept_mu<-accept_mu
  mcmc$accept_k<-accept_k
  mcmc$accept_theta_a<-accept_theta_a
  
  mcmc$nu_1<-mcmc$beta_v[1,1,]
  mcmc$nu_2<-mcmc$beta_v[1,2,]
  return(mcmc)
}

#' MCMC function for model with separate longitudinal and survival processes (SLS model)
#'
#' @param data_all A list of complete data. X0_inds is the baseline data for each patient, D is the dosage, Y is the longitudinal process, and Ts are the visit times,
#' id are the patient id's (ranging from 1 to the total number of patients), censor and surv_time are the censoring and survival times, and Npat is the total number of patients.
#' @param mcmc_settings A list for MCMC setup. burn.in is the number of total iterations to burn, ndisplay is number of iterations per which the display message will
#' appear, and Niter is the total number of iterations (including burn-in iterations).
#' @param seed Seed for running MCMC 
#'
#'
#'
#' @return 
#'   \item{beta_l}{Linear coefficient parameter in longitudinal submodel}
#'   \item{sigma2_l}{Error parameter in longitudinal submodel}
#'   \item{beta_d}{Linear coefficient parameter in dosing submodel}
#'   \item{sigma2_d}{Error parameter in dosing submodel}
#'   \item{nu_1}{Visitation intensity peak parameter in vistation submodel}
#'   \item{nu_2}{Visitation intensity shape parameter in vistation submodel}
#'   \item{mu}{Baseline visitation intensity parameter in vistation submodel}
#'   \item{beta_s1}{Creatinine-associated parameter in survival submodel}
#'   \item{beta_s2}{Dosage-associated parameter in survival submodel}
#'   \item{beta_s3}{Toxicity-associated parameter in survival submodel}
#'   \item{beta_s4}{Visitation-associated parameter in survival submodel}
#'   \item{h0}{Baseline hazard parameter in survival submodel}
#'   \item{shape}{Shape parameter in survival submodel}
#' 
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
#' mcmc_out_separate<-mcmc_separate(simulated_data,mcmc_settings,seed)
#' 
#' }
#' 
#' @export
#' 
mcmc_separate<-function(data_all,mcmc_settings,seed){
  Niter=mcmc_settings$Niter
  burn.in=mcmc_settings$burn.in
  ndisplay=mcmc_settings$ndisplay
  peak_dist=mcmc_settings$peak_dist
  D<-data_all$D
  Y<-data_all$Y
  Ts<-data_all$Ts
  id<-data_all$id
  censor<-data_all$censor
  surv_time<-data_all$surv_time
  Npat<-data_all$Npat
  
  #Prepping data
  inds_shrink_init<-shrink_data_inds(Npat,id,init=T)
  inds_shrink_end<-shrink_data_inds(Npat,id,init=F)
  inds_shrink_both<-shrink_data_inds_both(Npat,id)
  end_inds<-end_data_inds(Npat,id)
  inds_all<-inds_i(Npat,id)
  J<-obs_count(Npat,id)
  J_cumsum<-cumsum(J) 
  death_init_death<-data_init_death(Ts,surv_time,Npat,id)
  X0_l<-data_all$X0_inds[end_inds,]
  X0_d<-data_all$X0_inds[end_inds,]
  X0_inds<-data_all$X0_inds
  X0<-X<-X0_l
  Ncov_l<-ncol(X0_l)+4 
  Ncov_d<-ncol(X0_d)+2 
  Ncov_il<-3
  mcmc<-NULL 
  mcmc$Ncov_l=Ncov_l
  mcmc$Ncov_d=Ncov_d
  mcmc$Ncov_il=Ncov_il
  mcmc$beta_l<-array(NA, c(Ncov_l, Niter))
  mcmc$b_il<-array(NA, c(Ncov_il,Npat,Niter))
  mcmc$B<-array(NA,c(Ncov_il,Ncov_il,Niter))
  mcmc$beta_d<-array(NA, c(Ncov_d, Niter))
  mcmc$sigma2_l<-rep(NA, Niter)
  mcmc$sigma2_d<-rep(NA, Niter)
  mcmc$lambda_tox<-rep(NA, Niter)
  mcmc$k<-rep(NA,Niter)
  mcmc$beta_s<-rep(NA, Niter)
  mcmc$beta_sd<-rep(NA, Niter)
  mcmc$beta_sd_cum<-rep(NA, Niter)
  mcmc$beta_v<-array(NA, c(2,2, Niter))
  mcmc$theta_a<-array(NA, c(2,Niter))
  mcmc$mu<-rep(NA, Niter)
  mcmc$h0<-rep(NA, Niter)
  mcmc$shape<-rep(NA, Niter)
  mcmc$beta_alpha<-rep(NA, Niter)
  
  mcmc$beta_l[,1]<-init_beta_l(Ncov_l,0.5^2,D,X0_l,Y,Ts,id,Npat)
  mcmc$b_il[,,1]<-matrix(0,nrow=Ncov_il,ncol=Npat)
  b_il3<-mcmc$b_il[3,id,1]
  mcmc$lambda_tox[1]<-60
  mcmc$k[1]<-2
  mcmc$beta_s[1]<- 0
  mcmc$beta_sd[1]<- 0
  mcmc$beta_sd_cum[1]<- 0
  mcmc$h0[1]<- 0
  mcmc$shape[1]<- 1
  mcmc$beta_alpha[1]<- 0
  mcmc$beta_v[,,1]<-matrix(c(2,0,2,0),nrow=2)
  mcmc$beta_d[,1]<-0
  mcmc$beta_d[1,1]<-5
  mcmc$sigma2_l[1]<-0.5^2
  mcmc$sigma2_d[1]<-0.5^2
  mcmc$B[,,1]<-matrix(c(0.3^2,0,0,
                        0,0.2^2,0,
                        0,0,(5*10^(-4))^2),nrow=3)
  mcmc$theta_a[,1]<- c(9,-1.4)
  mcmc$mu[1]<- -4.5
  
  #hyperparameter values
  lambda1_l<-lambda2_l<-lambda1_d<-lambda2_d<-0.01
  a<-400
  b<-200
  c_B<-d_B<-0
  #MH settings + acceptance rates
  accept_beta_l<-0
  accept_beta_v<-matrix(c(0,0,0,0),nrow=2)
  accept_beta_s<-0 
  accept_beta_sd<-0
  accept_lambda_tox<-0
  accept_beta_sd_cum<-0 
  accept_mu<-0
  accept_theta_a<-c(0,0) 
  accept_h0<-0
  accept_shape<-0 
  accept_beta_alpha<-0
  accept_k<-0
  jump_v<-matrix(c(0.07,0,0.8,0),nrow=2)/3 
  jump_s<-0.06  
  jump_sd<-0.18
  jump_sd_cum<-0.18
  jump_lambda_tox<-7
  jump_mu<-0.05
  jump_theta_a<-c(0.4,0.08)/4.5
  jump_h0<-0.25 
  jump_k<-0.05 
  jump_shape<-0.04
  jump_beta_alpha<-0.4 
  miter=2  

  set.seed(seed) 
  for (iter in miter:Niter){
    if (iter%%ndisplay==0){
      cat("Iteration", iter,'of',Niter,'complete \n') 
      # print(mcmc$lambda_tox[iter-1])
      
    } 
    EY_out<-calculate_EY(mcmc$beta_l[,iter-1],mcmc$b_il[,,iter-1],D,Y,Ts,X0_l,Npat,inds_all)
    EY<-EY_out$EY
    mean_rand<-EY_out$rand
    int_all<-NA
    # beta_l<-mcmc$beta_l[,iter-1] ; beta_d<- mcmc$beta_d[,iter-1] ; beta_s<-mcmc$beta_s[iter-1] ; sigma2_l<-mcmc$sigma2_l[iter-1] ;
    # sigma2_d<-mcmc$sigma2_d[iter-1] ;  h0<- mcmc$h0[iter-1] ; b_il<- mcmc$b_il[,,iter-1] ; B<-mcmc$B[,,iter-1]
    # beta_s<-mcmc$beta_s[iter-1] ; h0<-mcmc$h0[iter-1]; shape<-mcmc$shape[iter-1];
    # theta_a<-mcmc$theta_a[,(iter-1)];beta_alpha<-mcmc$beta_alpha[iter-1];beta_sd=mcmc$beta_sd[iter-1];
    # beta_sd_cum=mcmc$beta_sd_cum[iter-1];Toxicity<-toxicity_calc_all(D,Ts,id,mcmc$lambda_tox[iter-1] ) ;lambda_tox=mcmc$lambda_tox[iter-1] ; k<-mcmc$k[iter-1]
    # k<-mcmc$k[iter-1]
    Toxicity<-toxicity_calc_all(D,Ts,id,mcmc$lambda_tox[iter-1],Npat )
    out<-update_beta_l_separate(mcmc$beta_l[,iter-1],mcmc$beta_d[,iter-1],mcmc$beta_s[iter-1],mcmc$sigma2_l[iter-1],
                       mcmc$sigma2_d[iter-1],D,X0_l,Y,Ts,id,Npat,int_all=NA,mcmc$h0[iter-1],mcmc$b_il[,,iter-1],
                       mean_rand,mcmc$B[,,iter-1],mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],
                       death_init_death,end_inds,inds_shrink_both, inds_shrink_init,inds_shrink_end,X0_inds,
                       mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],
                       Toxicity,mcmc$lambda_tox[iter-1],b_il3,inds_all,mcmc$k[iter-1],Ncov_l,surv_time,J_cumsum, J)
    mcmc$beta_l[,iter]<-out$beta_l
    mcmc$b_il[,,iter]<-out$b_il
    b_il3<-mcmc$b_il[3,id,iter]
    
    mcmc$B[,,iter]<-update_B(c_B,d_B,mcmc$b_il[,,iter],Npat)
    EY_out<-calculate_EY(mcmc$beta_l[,iter],mcmc$b_il[,,iter],D,Y,Ts,X0_l,Npat,inds_all)
    EY<-EY_out$EY
    mean_rand<-EY_out$rand
    mean_fixed<-EY_out$fixed
    # 
    # beta_s<-mcmc$beta_s[iter-1] ; beta_l<-mcmc$beta_l[,iter]; h0<-mcmc$h0[iter-1]; shape<-mcmc$shape[iter-1]; theta_a<-mcmc$theta_a[,(iter-1)];beta_alpha<-mcmc$beta_alpha[iter-1];beta_sd=mcmc$beta_sd[iter-1]; beta_sd_cum=mcmc$beta_sd_cum[iter-1];Toxicity<-toxicity_calc_all(D,Ts,id,mcmc$lambda_tox[iter-1] ) ;lambda_tox=mcmc$lambda_tox[iter-1]
    #  
    curr_loglik_surv=NA
    out_lambda_tox<-update_lambda_tox(mcmc$beta_s[iter-1],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,
                                      jump_lambda_tox,mcmc$h0[iter-1],mean_rand,Y,mcmc$shape[iter-1],
                                      censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,
                                      inds_shrink_end,mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],Toxicity,mcmc$lambda_tox[iter-1],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$lambda_tox[iter]<-out_lambda_tox$lambda_tox
    jump_lambda_tox<- sqrt(out_lambda_tox$lambda_tox)*5
    Toxicity=out_lambda_tox$Toxicity
    curr_loglik_surv<-out_lambda_tox$curr_loglik_surv
    
    out_s<-update_beta_s(mcmc$beta_s[iter-1],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_s,mcmc$h0[iter-1],mean_rand,Y,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_s[iter]<-out_s$beta_s
    curr_loglik_surv<-out_s$curr_loglik_surv
    
    out_sd<-update_beta_sd(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_sd,mcmc$h0[iter-1],mean_rand,Y,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter-1], mcmc$beta_sd_cum[iter-1],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_sd[iter]<-out_sd$beta_sd
    curr_loglik_surv<-out_sd$curr_loglik_surv
    
    out_sd_cum<-update_beta_sd_cum(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_sd_cum,mcmc$h0[iter-1],mean_rand,Y,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter-1],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_sd_cum[iter]<-out_sd_cum$beta_sd_cum
    curr_loglik_surv<-out_sd_cum$curr_loglik_surv
    
    out_h0<-update_h0(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_h0,mcmc$h0[iter-1],mean_rand,Y,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$h0[iter]<-out_h0$h0
    curr_loglik_surv<-out_h0$curr_loglik_surv
    
    out_shape<-update_shape(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_shape,mcmc$h0[iter],mean_rand,Y,mcmc$shape[iter-1],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity ,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$shape[iter]<-out_shape$shape
    curr_loglik_surv<-out_shape$curr_loglik_surv
    
    out_beta_alpha<-update_beta_alpha(mcmc$beta_s[iter],mcmc$beta_l[,iter],D,Y,Ts,id,Npat,surv_time,jump_beta_alpha,mcmc$h0[iter],mean_rand,Y,mcmc$shape[iter],censor,mcmc$theta_a[,(iter-1)],mcmc$beta_alpha[iter-1],death_init_death,end_inds,inds_shrink_end,mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity,mcmc$lambda_tox[iter],b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l)
    mcmc$beta_alpha[iter]<-out_beta_alpha$beta_alpha
    curr_loglik_surv<-out_beta_alpha$curr_loglik_surv
    
    out_v<-update_beta_v(mcmc$beta_v[,,iter-1],Y,Ts,mcmc$mu[iter-1],mcmc$theta_a[,iter-1],Npat,jump_v,id,inds_shrink_init,inds_shrink_end,mcmc$k[iter-1],peak_dist)
    mcmc$beta_v[,,iter]<-out_v$beta_v
    
    out_mu<-update_mu(mcmc$beta_v[,,iter],Y,Ts,mcmc$mu[iter-1],mcmc$theta_a[,iter-1],Npat,jump_mu,id,inds_shrink_init,inds_shrink_end,mcmc$k[iter-1],peak_dist)
    mcmc$mu[iter]<-out_mu$mu
    
    out_theta_a<-update_theta_a(mcmc$beta_v[,,iter],Y,Ts,mcmc$mu[iter],mcmc$theta_a[,iter-1],Npat,jump_theta_a,mcmc$beta_alpha[iter],mcmc$beta_s[iter],
                                surv_time,mcmc$beta_l[,iter],D,id,mcmc$h0[iter],mean_rand,Y,mcmc$shape[iter],censor,inds_shrink_init,inds_shrink_end,death_init_death,end_inds,
                                mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity,mcmc$lambda_tox[iter],inds_shrink_both,b_il3,curr_loglik_surv,mcmc$k[iter-1],Ncov_l,peak_dist)
    mcmc$theta_a[,iter]<-out_theta_a$theta_a
    curr_loglik_surv<-out_theta_a$curr_loglik_surv
    
    out_k<-update_k(mcmc$beta_v[,,iter],Y,Ts,mcmc$mu[iter],mcmc$theta_a[,iter],Npat,jump_theta_a,mcmc$beta_alpha[iter],mcmc$beta_s[iter],
                    surv_time,mcmc$beta_l[,iter],D,id,mcmc$h0[iter],mean_rand,Y,mcmc$shape[iter],censor,inds_shrink_init,inds_shrink_end,death_init_death,end_inds,
                    mcmc$beta_sd[iter], mcmc$beta_sd_cum[iter],Toxicity,mcmc$lambda_tox[iter],inds_shrink_both,b_il3,curr_loglik_surv,mcmc$k[iter-1],jump_k,a,b,Ncov_l,peak_dist)
    mcmc$k[iter]<-out_k$k
    
    mcmc$beta_d[,iter]<-update_beta_d(mcmc$beta_d[,iter-1],mcmc$sigma2_d[iter-1],X0_inds,Y,D,Ts,id,int_all,Npat)
    
    mcmc$sigma2_l[iter]<-update_sigma2_l(lambda1_l,lambda2_l,EY,Y,Npat,id,mcmc$B[iter],inds_shrink_init)
    mcmc$sigma2_d[iter]<-update_sigma2_d(lambda1_d,lambda2_d,Y,Npat,id,int_all,mcmc$beta_d[,iter],mcmc$B[iter],X0_inds,D)
    
    if (iter>burn.in){
      accept_beta_l<-accept_beta_l+out$accept_fixed
      accept_beta_v<-accept_beta_v+out_v$accept 
      accept_beta_s<-accept_beta_s+out_s$accept
      accept_lambda_tox<-accept_lambda_tox+out_lambda_tox$accept
      
      accept_beta_sd<-accept_beta_sd+out_sd$accept 
      accept_beta_sd_cum<-accept_beta_sd_cum+out_sd_cum$accept
      accept_h0<-accept_h0+out_h0$accept
      accept_shape<-accept_shape+out_shape$accept
      accept_beta_alpha<-accept_beta_alpha+out_beta_alpha$accept
      accept_mu<-accept_mu+out_mu$accept
      accept_k<-accept_k+out_k$accept
      accept_theta_a<-accept_theta_a+out_theta_a$accept
    }
    
    miter<-iter
  }
  mcmc$accept_beta_l<-accept_beta_l
  mcmc$accept_beta_v<-accept_beta_v
  mcmc$accept_beta_s<-accept_beta_s
  mcmc$accept_lambda_tox<-accept_lambda_tox
  mcmc$accept_beta_sd<-accept_beta_sd
  mcmc$accept_beta_sd_cum<-accept_beta_sd_cum
  mcmc$accept_h0<-accept_h0
  mcmc$accept_shape<-accept_shape
  mcmc$accept_beta_alpha<-accept_beta_alpha
  mcmc$accept_mu<-accept_mu
  mcmc$accept_k<-accept_k
  mcmc$accept_theta_a<-accept_theta_a 
  return(mcmc)
}

