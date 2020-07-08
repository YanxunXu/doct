#' @import MASS
#' @import pracma
#' @import mvtnorm 
#' @importFrom LaplacesDemon rinvwishart
#' @import stats
#' @import GoFKernel
#' @useDynLib doct
#' @importFrom Rcpp evalCpp
#' @importFrom survival survfit


int_func<-function(ti,t_r,y_r,beta_v,theta_a,mu,k){ 
  nu=exp(as.vector(t(beta_v[,1])%*%c(1,y_r)))
  kappa=exp(as.vector(t(beta_v[,2])%*%c(1,y_r)))+1
  alpha=k/(1+exp(as.numeric(theta_a%*%t(cbind(1,y_r)))))
  
  return(exp(mu)+alpha*dgamma((ti-t_r),shape=kappa,rate=((kappa-1)/nu) ))
} 



random_visit_times_Rcpp = function(int_func,t_r,y_r,beta_v,theta_a,mu,k,max_time=10000){
  alpha=sigmoid_k(as.numeric((theta_a)%*%c(1,y_r)),k)  #Use this instead of exp(-t_r/theta_a) for now.
  cumulative_density_fn = function(t) 1 - exp(-integrate_intensity(mu, beta_v[1,1], beta_v[1,2],alpha,0,t))
  inverse_cumulative_density_fn = inverse(cumulative_density_fn, 0, max_time)
  u<-runif(1)
  t_new<-max(inverse_cumulative_density_fn(u),0.001)+t_r
  out<-NULL
  out$t_new<-t_new
  out$int<-int_func(t_new,t_r,y_r,beta_v,theta_a,mu,k)
  return(out)
}



haz_func_mixed_recent<-function(ti,t_r,beta_s,h0,e_y,shape,beta_alpha,alpha,beta_sd,beta_sd_cum,di,tox,lambda_tox,beta_l,b_il3){
  tox_updated<- ((tox*exp(-(ti-t_r)/lambda_tox)) + (1-exp(-(ti-t_r)/lambda_tox ) ) * di )
  Ncov<-length(beta_l)
  e_y_updated<-e_y+(ti-t_r)*(beta_l[Ncov-1]+b_il3)+(ti-t_r)^2*beta_l[Ncov]
  haz<-(shape)*exp(-(beta_s*e_y_updated+beta_alpha*alpha+beta_sd*di+beta_sd_cum*tox_updated+h0))*ti^(shape-1)
  haz[is.infinite(haz)]=0
  return(haz)
}


integrate_from_0_surv_mixed_recent = function(fn, t,t_r,beta_s,h0,e_y,shape,beta_alpha,alpha,beta_sd,beta_sd_cum,di,tox,lambda_tox,beta_l,b_il3=b_il3){
  int_fn = function(t) integrate(fn, t_r, t,t_r=t_r,beta_s=beta_s,h0=h0,e_y=e_y,shape=shape,beta_alpha=beta_alpha,alpha=alpha,
                                 beta_sd=beta_sd,beta_sd_cum=beta_sd_cum,di=di,tox=tox,lambda_tox=lambda_tox,beta_l=beta_l,b_il3=b_il3)
  result = sapply(t, int_fn)
  value  = unlist(result["value",]) 
  msg    = unlist(result["message",])
  value[which(msg != "OK")] = NA
  return(value)
}
 
random_death_times_mixed_recent = function(haz_func_mixed_recent,t_max,t_r,beta_s,h0,e_y,shape,beta_alpha,alpha,beta_sd,beta_sd_cum,di,tox,lambda_tox,beta_l,b_il3){
  death<-NULL
  cumulative_density_fn_surv_mixed_recent = function(t) 1 - exp(-integrate_from_0_surv_mixed_recent(haz_func_mixed_recent, t,t_r,beta_s,h0,e_y,shape,beta_alpha,alpha,beta_sd,beta_sd_cum,di,tox,lambda_tox,beta_l,b_il3))
  inverse_cumulative_density_fn_surv_mixed_recent = inverse(cumulative_density_fn_surv_mixed_recent, t_r, t_max)
  u<-runif(1)
  if (cumulative_density_fn_surv_mixed_recent(t_max)<u){
    death$ind<-0
  }
  else{
    t_new<-inverse_cumulative_density_fn_surv_mixed_recent(u)
    if (t_new<t_max){
      death$ind<-1
      death$time<-t_new
    } else{
      death$ind<-0
    }
  }
  return(death)
}

toxicity_calc<-function(d_i,t_i,t_r,lambda_tox,tox_recent){
  tox<-tox_recent*exp(-(t_i-t_r)/lambda_tox)
  weight=(1-exp(-(t_i-t_r)/lambda_tox ) )
  tox=tox+d_i*weight
  return(tox)
}


toxicity_calc_i<-function(D,Ts,lambda_tox){
  Toxicity=numeric(length(D))
  d_i<-D
  t_i<-Ts
  n<-length(d_i)
  tox_tmp=numeric(n)
  for (len in 2:n){
    tox<-tox_tmp[len-1]*exp(-(t_i[len]-t_i[len-1])/lambda_tox)
    weight=(1-exp(-(t_i[len]-t_i[len-1])/lambda_tox ) )
    tox=tox+d_i[len-1]*weight
    tox_tmp[len]<-tox
  }
  return(tox_tmp)
}



toxicity_calc_all<-function(D,Ts,id,lambda_tox,Npat){
  Toxicity=numeric(length(D))
  for (i in 1:Npat){
    inds<-which(id==i)
    len<-length(inds)
    d_i<-D[inds]
    t_i<-Ts[inds]
    n<-length(d_i)
    tox_tmp=numeric(n)
    for (len in 2:n){
      tox<-tox_tmp[len-1]*exp(-(t_i[len]-t_i[len-1])/lambda_tox)
      weight=(1-exp(-(t_i[len]-t_i[len-1])/lambda_tox ) )
      tox=tox+d_i[len-1]*weight
      tox_tmp[len]<-tox
    }
    Toxicity[inds]=tox_tmp
  } 

  return(Toxicity)
}

#' Generate Simulated Data
#'
#' @param Npat Number of patients to simulate.
#' @param beta_l Linear coefficient parameter in longitudinal submodel
#' @param sigma2_l Error parameter in longitudinal submodel
#' @param beta_d Linear coefficient parameter in dosing submodel
#' @param sigma2_d Error parameter in dosing submodel
#' @param nu_1 Visitation intensity peak parameter in vistation submodel
#' @param nu_2 Visitation intensity shape parameter in vistation submodel
#' @param mu Baseline visitation intensity parameter in vistation submodel
#' @param beta_s1 Creatinine-associated parameter in survival submodel
#' @param beta_s2 Dosage-associated parameter in survival submodel
#' @param beta_s3 Toxicity-associated parameter in survival submodel
#' @param beta_s4 Visitation-associated parameter in survival submodel
#' @param h0 Baseline hazard parameter in survival submodel
#' @param s Shape parameter in survival submodel
#'
#' @return List of data. X0_inds is the baseline data for each patient, D is the dosage, Y is the longitudinal process, and Ts are the visit times,
#' id are the patient id's (ranging from 1 to the total number of patients), censor and surv_time are the censoring and survival times, and Npat is the total number of patients.
#'
#'@examples
#' \dontrun{
#' #Simulate Data
#' seed=203
#' set.seed(seed) 
#' data_all<-Generate_Simulated_Data(10)
#' }
#' @export
#' 

Generate_Simulated_Data<-function(Npat,beta_l=c(5.3,0.1,0.3,0.4, 0.25, -1*10^(-4),3*10^(-8)),sigma2_l=0.1^2,
                        beta_d=c(1,0.2,0.15,0.2,0.15),
                        sigma2_d=0.3^2,nu_1=2.5,nu_2=1.5,mu=-4.8,beta_s1=1,
                        beta_s2=0.9,beta_s3=-0.75,beta_s4=-5,h0=5,s=1.05){ 
  #Create data: tacrotl, creat, survival
  
  # beta_l<-c(5.3,0.1,0.3,0.4, 0.25, -1*10^(-4),3*10^(-8))
  # sigma2_l<-0.1^2
  # 
  # beta_d<-c(1,0.2,0.15,0.2,0.15)
  # sigma2_d<-0.3^2 

  B<-matrix(c(0.2^2,0,0,
              0,0.07^2, 0,
              0,0,(1*10^(-4))^2),nrow=3)
  beta_v<-matrix(c(nu_1,0,nu_2,0.00),nrow=2) 
  
  theta_a<-c(9.5,-1.5)
  # mu<- -4.8
  
  beta_s=beta_s1 #prev
  beta_sd=beta_s2
  beta_sd_cum=beta_s3
  # h0=5
  shape<-s
  lambda_tox<-50
  beta_alpha=beta_s4
  k<-2
  tot_censor<-0
  for (i in 1:(Npat))
  {
    t_i<-0
    DGF_i<-rbinom(1,1,0.4)
    ageD_i<-rnorm(1,0,1)
    diab_i<-rbinom(1,1,0.2)
    type_i<-rbinom(1,1,0.8)
    BMI_i<-rnorm(1,0,1)
    death_ind=0
    mean_y_init<-5
    y_i<-rnorm(1,mean_y_init,sqrt(sigma2_l))
    
    cov_d<-c(1,y_i,ageD_i,DGF_i, BMI_i)
    d_i<-rnorm(1,t(beta_d)%*%cov_d,sqrt(sigma2_d))
    b_il<-mvrnorm(1,mu=c(0,0,0),Sigma=B)
    tox_i<-0
    mean_rand<-as.numeric(b_il[1])
    mean_rand_i<-mean_rand
    mean_fixed_i<-mean_y_init
    ey_i<-mean_y_init
    ints<-NA
    j=1
    alphas<-c()
    censor_time<-rweibull(1,3,8000)
    censor=1
    death_i<-censor_time
    while(death_ind==0 && censor==1){
      j<-j+1
      hawk_out<-random_visit_times_Rcpp(int_func,t_i[j-1],y_i[j-1],beta_v,theta_a,mu,k)
      t_ij<-hawk_out$t_new 
      if(t_ij>censor_time){
        censor=0
        t_i_all<-t_i
        tot_censor=tot_censor+1
      } else{
        t_i<-c(t_i,t_ij)
        alpha=sigmoid_k(as.vector(t(theta_a)%*%c(1,y_i[j-1])),k)
        
        alphas<-c(alphas,alpha)
        death<-random_death_times_mixed_recent(haz_func_mixed_recent,t_ij,t_i[j-1],beta_s,h0,ey_i[j-1],shape,beta_alpha,alpha,beta_sd,beta_sd_cum,d_i[j-1],tox_i[j-1],lambda_tox,beta_l,b_il[3])

        if(death$ind==1 & j>2){
          death_ind=1 
          t_i_all<-t_i
          t_i<-t_i[-j]
          death_i<-death$time
        } else{
          ints<-c(ints,hawk_out$int)
          # x_ij<-rnorm(1,0.5,1)
          # x_ij<-1
          # x_i<-c(x_i,x_ij)
          w_bmi<-0.1
          cov_l<-c(1,d_i[j-1],ageD_i,DGF_i, BMI_i, t_i[j],t_i[j]^2)
          cov_il<-c(1,d_i[j-1],t_i[j])
          mean_fixed<-t(beta_l)%*%cov_l
          mean_rand<-as.numeric(t(b_il)%*%cov_il)

          e_ij<-mean_rand+mean_fixed
          mean_rand_i<-c(mean_rand_i,mean_rand)
          mean_fixed_i<-c(mean_fixed_i,mean_fixed)
          
          ey_i<-c(ey_i,e_ij)
          y_ij<-rnorm(1,e_ij,sqrt(sigma2_l))
          y_i<-c(y_i,y_ij)
          cov_d<-c(1,y_ij,ageD_i,DGF_i, BMI_i)
          
          
          d_ij<-rnorm(1,t(beta_d)%*%cov_d,sqrt(sigma2_d))
          d_i<-c(d_i,d_ij)
          tox_i<-c(tox_i,toxicity_calc(d_i[j-1],t_i[j],t_i[j-1],lambda_tox,tox_i[j-1]))
        }
      }
    }
    # plot(t_i,y_i,ty='l')
    # plot(t_i,d_i,ty='l')
    Data_Reg_temp<-NULL
    Data_Reg_temp$ageD<-rep(ageD_i,length(t_i))
    Data_Reg_temp$Y<-y_i
    Data_Reg_temp$D<-d_i
    Data_Reg_temp$Toxicity<-tox_i
    Data_Reg_temp$DGF<-rep(DGF_i,length(t_i))
    # Data_Reg_temp$diab<-rep(diab_i,length(t_i))
    # Data_Reg_temp$typedonor<-rep(type_i,length(t_i))
    Data_Reg_temp$BMI<-rep(BMI_i,length(t_i))
    
    
    Data_Reg_temp$Ts<-t_i
    Data_Reg_temp$Ts_death<-t_i_all
    
    Data_Reg_temp$id<-rep(i,length(t_i))
    Data_Reg_temp$id_t<-rep(i,length(t_i_all))
    
    # Data_Reg_temp$ints<-ints
    Data_Reg_temp$EY<-ey_i
    Data_Reg_temp$mean_fixed<-mean_fixed_i
    Data_Reg_temp$mean_rand<-mean_rand_i
    
    Data_Reg_temp$J<-j-1
    Data_Reg_temp$surv_time<-rep(death_i,length(t_i))
    Data_Reg_temp$surv_time_i<-death_i
    
    Data_Reg_temp$censor<-rep(censor,length(t_i))
    # Data_Reg_temp$censor_time<-rep(censor_time,length(t_i))
    
    Data_Reg_temp$b_il<-b_il
    # Data_Reg_temp$alphas<-alphas
    
    if (i==1){
      Data_Reg<-Data_Reg_temp
    }else{
      Data_Reg<-mapply(c, Data_Reg, Data_Reg_temp, SIMPLIFY=FALSE)
    }
  }
  
  Data_Reg$Npat<-Npat
  Data_Reg$total_censored<-tot_censor
  Data_Reg$b_il<-matrix(Data_Reg$b_il,nrow=3)
  Data_Reg$beta_l=beta_l
  Data_Reg$B=B
  Data_Reg$sigma2_l=sigma2_l
  Data_Reg$beta_d=beta_d
  Data_Reg$sigma2_d=sigma2_d
  Data_Reg$beta_v=beta_v
  Data_Reg$theta_a=theta_a
  Data_Reg$mu=mu
  Data_Reg$beta_s=beta_s
  Data_Reg$beta_sd=beta_sd
  Data_Reg$beta_sd_cum=beta_sd_cum
  Data_Reg$h0=h0
  Data_Reg$beta_alpha=beta_alpha
  Data_Reg$shape=shape
  Data_Reg$lambda_tox=lambda_tox
  Data_Reg$k=k
  Data_Reg$X0_inds<-cbind(Data_Reg$ageD,Data_Reg$DGF,Data_Reg$BMI)
  
  return(Data_Reg)
}


#' Simulated dataset for 500 patients.
#'
#'
#' @format List of data. X0_inds is the baseline data for each patient, D is the dosage, Y is the longitudinal process, and Ts are the visit times,
#' id are the patient id's (ranging from 1 to the total number of patients), censor and surv_time are the censoring and survival times, and Npat is the total number of patients.
"simulated_data"

