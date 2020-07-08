// [[Rcpp::depends(RcppNumerical)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]] 

#include <RcppArmadillo.h>
#include <RcppNumerical.h>
#include <RcppEigen.h>

using namespace Rcpp; 
using namespace Numer;
using namespace arma;

class intensity: public Numer::Func {
private:
  double mu;
  double beta_v11; 
  double beta_v12;
  double alpha;
public:
  intensity(double mu_, double beta_v11_, double beta_v12_, double alpha_) : mu(mu_), beta_v11(beta_v11_), beta_v12(beta_v12_), alpha(alpha_) {}
  double operator()(const double& x) const {
    double nu=std::exp(beta_v11) ;
    double kappa=std::exp(beta_v12)+1 ;
    return (std::exp(mu) + alpha*pow(x,(kappa-1))*std::exp(-(kappa-1)/nu*(x))*pow( ((kappa-1)/nu),(kappa) )/std::tgamma(kappa) );
    // return (std::exp(mu) + alpha*R::dgamma(x,kappa,nu/(kappa-1),0));
    // return (std::exp(mu) + alpha*Rcpp::dgamma(x,kappa,nu/(kappa-1),false));
    
  }
}; 


class haz_func: public Numer::Func {
private:
  double shape;
  double beta_s; 
  double e_y;
  double alpha;
  double beta_alpha;
  double h0;
  double beta_sd;
  double beta_sd_cum;
  double di;
  double tox;
  double lambda_tox;
  double t_r;
  arma::vec beta_l;
  double b_il3;
  int Ncov_l;
public:
  haz_func(double shape_, double beta_s_, double e_y_, double alpha_, 
           double beta_alpha_, double h0_, double beta_sd_,
           double beta_sd_cum_, double di_, double tox_, double lambda_tox_, double t_r_, arma::vec beta_l_, double b_il3_, int Ncov_l_): shape(shape_), 
           beta_s(beta_s_), e_y(e_y_), alpha(alpha_), 
           beta_alpha(beta_alpha_), h0(h0_), beta_sd(beta_sd_), 
           beta_sd_cum(beta_sd_cum_), di(di_), tox(tox_), lambda_tox(lambda_tox_), t_r(t_r_) , beta_l(beta_l_), b_il3(b_il3_), Ncov_l(Ncov_l_){}
  
  double operator()(const double& x) const {
    // return ((shape)*exp(-(beta_s*e_y+beta_alpha*alpha+h0))*(pow(x,(shape-1))) );
    double e_y_updated=e_y+(x-t_r)*(beta_l[(Ncov_l-2)]+b_il3)+ pow((x-t_r),2)*beta_l[(Ncov_l-1)];
    
    // double tox_updated= ((100*tox*exp(-(x-t_r)/lambda_tox)) + (1-exp(-(x-t_r)/lambda_tox ) )*lambda_tox * di )/100;
    double tox_updated= ((tox*exp(-(x-t_r)/lambda_tox)) + (1-exp(-(x-t_r)/lambda_tox ) ) * di );
    
    // double tox_updated= tox;
    
    return ((shape)*exp(-(beta_s*e_y_updated + beta_alpha*alpha + beta_sd*di + beta_sd_cum*tox_updated + h0))*(pow(x,(shape-1)))  );
    
    
  } 
}; 

// [[Rcpp::export]]
double haz_func_gen(double shape, double beta_s, double e_y, double alpha, 
                    double beta_alpha, double h0, double beta_sd,
                    double beta_sd_cum, double di, double tox, double lambda_tox, double x, double t_r, arma::vec beta_l, double b_il3, int Ncov_l) {
  
  double e_y_updated=e_y+(x-t_r)*(beta_l[(Ncov_l-2)]+b_il3)+ pow((x-t_r),2)*beta_l[(Ncov_l-1)];
  double tox_updated= ((tox*exp(-(x-t_r)/lambda_tox)) + (1-exp(-(x-t_r)/lambda_tox ) ) * di );
  return ((shape)*exp(-(beta_s*e_y_updated + beta_alpha*alpha + beta_sd*di + beta_sd_cum*tox_updated + h0))*(pow(x,(shape-1)))  );
} 




// [[Rcpp::export]]
double test_haz(double shape, double beta_s, double e_y, double alpha, double beta_alpha, double h0, double ti,double beta_sd,
                double beta_sd_cum, double di, double tox, double lambda_tox, double t_r) {
  double tox_updated= ((100*tox*exp(-(ti-t_r)/lambda_tox)) + (1-exp(-(ti-t_r)/lambda_tox ) )*lambda_tox * di )/100;
  return ((shape)*exp(-(beta_s*e_y + beta_alpha*alpha + beta_sd*di + beta_sd_cum*tox_updated + h0))*(pow(ti,(shape-1))) );
} 

// [[Rcpp::export]]
double test_intensity(double x, double mu,double beta_v11, double beta_v12, double alpha) {
  double nu=std::exp(beta_v11) ; 
  double kappa=std::exp(beta_v12)+1 ;
  
  return (std::exp(mu) + alpha*R::dgamma(x,kappa,nu/(kappa-1),0));
  // return (std::exp(mu) + alpha*Rcpp::dgamma(x,kappa,nu/(kappa-1),0) );
  
} 


// [[Rcpp::export]]
double integrate_intensity(const double &mu,const double &beta_v11,const double &beta_v12,const double &alpha,
                                        const double &lower, const double &upper) {
  intensity function(mu,beta_v11,beta_v12,alpha);
  double err_est;
  int err_code;
  int subdiv=25;
  double eps_abs= pow(10,-1);
  double eps_rel= pow(10,-1);
  
  const double result = Numer::integrate(function, lower, upper, err_est, err_code, subdiv,
                                         eps_abs,eps_rel, Integrator<double>::GaussKronrod15);
  return  result;
  // ,
  // Rcpp::Named("error") = err_est
}

// [[Rcpp::export]] 
double integrate_hazard(const double &shape,const double &beta_s,const double &e_y,const double &alpha,const double &beta_alpha,const double &h0,
                                     const double &beta_sd,const double &beta_sd_cum, const double &di, const double &tox,const double &lambda_tox,const double &t_r,
                                     arma::vec beta_l, const double b_il3 , const double &lower, const double &upper, int Ncov_l) {
  haz_func function(shape,beta_s,e_y,alpha,beta_alpha,h0,beta_sd,beta_sd_cum,di,tox,lambda_tox,t_r,beta_l,b_il3,Ncov_l);
  double err_est;
  int err_code;
  int subdiv=25;
  double eps_abs= pow(10,-1);
  double eps_rel= pow(10,-1);
  const double result = Numer::integrate(function, lower, upper, err_est, err_code, subdiv,
                                         eps_abs,eps_rel, Integrator<double>::GaussKronrod15);
  // const double result = Numer::integrate(function, lower, upper, err_est, err_code);
  // const double result = Numer::integrate(function, lower, upper, err_est, err_code,subdiv,eps_abs,eps_rel);
  
  return  result;
  // ,
  // Rcpp::Named("error") = err_est
}


// [[Rcpp::export]]
double loglike_surv_Rcpp(const double &shape,const double &beta_s,const double &beta_alpha,const double &h0,
                         const double &beta_sd,const double &beta_sd_cum, arma::vec D,const double &lambda_tox,
                         arma::vec inds_shrink_end, arma::vec end_inds,arma::vec EY,
                         arma::vec alphas,arma::vec Ts, arma::vec Toxicity,arma::vec death_time, 
                         int Npat, int Nobs, arma::vec beta_l, arma::vec b_il3, int Ncov_l) {
  double log_surv_all=0;
  for (int i=0; i< Npat; ++i){
    int j=end_inds[i]-1;
    // std::cout << j << std::endl; 
    
    log_surv_all=log_surv_all- integrate_hazard(shape,beta_s,EY[j],alphas[j],beta_alpha,h0,beta_sd,beta_sd_cum, D[j], Toxicity[j],lambda_tox,Ts[j],beta_l,b_il3[j],
                                                           Ts[j], death_time[j],Ncov_l) ;
  }
  for (int i=0; i< Nobs ; ++i){
    int j=inds_shrink_end[i]-1;
    log_surv_all=log_surv_all- integrate_hazard(shape,beta_s,EY[j],alphas[j],beta_alpha,h0,beta_sd,beta_sd_cum, D[j], Toxicity[j],lambda_tox,Ts[j],beta_l,b_il3[j],
                                                           Ts[j], Ts[j+1],Ncov_l) ;
  }
  return(log_surv_all);
}


// [[Rcpp::export]]
double loglike_surv_i_Rcpp(const double &shape,const double &beta_s,const double &beta_alpha,const double &h0,
                           const double &beta_sd,const double &beta_sd_cum, arma::vec D,const double &lambda_tox, arma::vec EY,
                           arma::vec alphas,arma::vec Ts, arma::vec Toxicity, double death_time, int Npat, int Nobs, arma::vec beta_l, double b_il3, int Ncov_l) {
  double log_surv_all=0;
  int j=Nobs-1;
  log_surv_all=log_surv_all- integrate_hazard(shape,beta_s,EY(j),alphas(j),beta_alpha,h0,beta_sd,beta_sd_cum, D(j), Toxicity(j),lambda_tox,Ts(j),beta_l,b_il3,
                                                         Ts(j), death_time,Ncov_l) ;
  // std::cout << log_surv_all << std::endl;
  // for (int j=1; j< (Nobs-1) ; ++j){
  for (int j=0; j< (Nobs-1) ; ++j){
    log_surv_all=log_surv_all- integrate_hazard(shape,beta_s,EY(j),alphas(j),beta_alpha,h0,beta_sd,beta_sd_cum, D(j), Toxicity(j),lambda_tox,Ts(j),beta_l,b_il3,
                                                           Ts(j), Ts(j+1),Ncov_l) ;
    // std::cout << log_surv_all << std::endl;
  }
  return(log_surv_all);
}

// [[Rcpp::export]]
arma::vec test(arma::mat y, arma::vec theta_a) {
  arma::vec alphas=1/(1+exp(y*theta_a));
  return(alphas);
  }
// [[Rcpp::export]]
List update_b_il(arma::vec J_cumsum, arma::vec J, int Npat, arma::vec D, arma::vec Y, arma::vec Ts, 
                 arma::vec curr_mean_fixed, arma::vec censor, arma::vec Toxicity, arma::mat b_il, arma::vec death_time,
                 arma::mat B, double sigma2_l,arma::vec curr_EY, arma::vec theta_a, double beta_s,
                 double shape, double h0, double beta_alpha, double beta_sd, double beta_sd_cum, double lambda_tox, arma::vec beta_l, double k, int Ncov_l) {
  int min_it=0;
  arma::vec loglike_allpat(Npat);
  arma::vec accept_rand=arma::zeros(Npat);
  for (int i=0; i<Npat ; ++i){
    int max_it=J_cumsum[i]-1;
    int j_i=J(i);
    arma::vec d_i=D.subvec(min_it,max_it);
    arma::vec y_i=Y.subvec(min_it,max_it);
    arma::vec b_il_i=b_il.col(i);
    arma::vec tox_i=Toxicity.subvec(min_it,max_it);
    double censor_stat=censor(min_it);
    double death_i=death_time(min_it);
    arma::vec mean_fixed_i=curr_mean_fixed.subvec(min_it,max_it);
    arma::vec time=Ts.subvec(min_it,max_it);
    arma::vec tmp=arma::ones(j_i-1,1);
    arma::mat cov_il= arma::join_rows(tmp,d_i.subvec(0,(j_i-2)), time.subvec(1,(j_i-1)));
    arma::mat Sigma=inv(inv(B)+(cov_il.t()*cov_il)/sigma2_l);
    arma::vec mu=1/sigma2_l*Sigma*(cov_il.t())*(y_i.subvec(1,j_i-1)-mean_fixed_i.subvec(1,j_i-1));

    arma::vec prop_b_il=mvnrnd(mu,Sigma);
    arma::vec ey_i_curr=curr_EY.subvec(min_it,max_it);
    arma::vec mean_rand_prop_i=arma::zeros(j_i);
    mean_rand_prop_i.subvec(1,(j_i-1))=cov_il*prop_b_il;
    arma::vec ey_i_prop= mean_fixed_i+mean_rand_prop_i;
    
    arma::vec tmp1=arma::ones(j_i,1);
    arma::mat y_i1= arma::join_rows(tmp1,y_i);
    arma::vec alphas=k/(1+exp(y_i1*theta_a));
    // alphas.print();
    // std::cout <<ey_i_prop[(j_i-1)] << std::endl;
    // std::cout <<ey_i_curr[(j_i-1)] << std::endl;
    
    double loghaz_all=log(haz_func_gen(shape, beta_s,ey_i_prop[(j_i-1)], alphas[(j_i-1)], beta_alpha, h0, beta_sd,
                                       beta_sd_cum,d_i[(j_i-1)],tox_i[(j_i-1)],lambda_tox, death_i,time[(j_i-1)],  beta_l,  b_il_i[3], Ncov_l))-
                                         log(haz_func_gen(shape, beta_s,ey_i_curr[(j_i-1)], alphas[(j_i-1)], beta_alpha, h0, beta_sd,
                                                          beta_sd_cum,d_i[(j_i-1)],tox_i[(j_i-1)],lambda_tox, death_i,time[(j_i-1)],  beta_l,  b_il_i[3], Ncov_l));
    // double loghaz_all=sum(log(haz_func(death_i,time[(j_i-1)],beta_s,h0,ey_i_prop[(j_i-1)],shape,beta_alpha,alphas[(j_i-1)],
    //                                           beta_sd,beta_sd_cum,d_i[(j_i-1)],tox_i[(j_i-1)],lambda_tox, beta_l,b_il_i[3]))*censor_stat);
    
    double loglike=loglike_surv_i_Rcpp(shape, beta_s, beta_alpha,h0,
                                       beta_sd, beta_sd_cum, d_i,lambda_tox,ey_i_prop,
                                       alphas, time,  tox_i, death_i, Npat,  j_i,  beta_l,  b_il_i[3],Ncov_l)-
                                         loglike_surv_i_Rcpp(shape, beta_s, beta_alpha,h0,
                                                             beta_sd, beta_sd_cum, d_i,lambda_tox,ey_i_curr,
                                                             alphas, time,  tox_i, death_i, Npat,  j_i,  beta_l,  b_il_i[3],Ncov_l);
    double loglike_i=loghaz_all*censor_stat+loglike;
    loglike_allpat[i]=loglike_i;
    if(as<double>(Rcpp::runif(1)) <exp(loglike_i)){
      b_il.col(i)=prop_b_il;
      accept_rand(i)=accept_rand[i]+1;
    }
    
    // tmp.insert_cols(1, d_i.subvec(0,(j_i-2)));
    // cov_il.print();
    
    min_it=max_it+1;
  }
  List output;
  output["b_il"]=b_il;
  output["accept_rand"]=accept_rand;
  
  return(output);
}


// [[Rcpp::export]]
double CDF_intensity(double mu,double beta_v11, double beta_v12, double alpha,
                     double lower,  double upper) {
  intensity function(mu,beta_v11,beta_v12,alpha);
  double err_est;
  int err_code;
  const double result = 1-std::exp(-Numer::integrate(function, lower, upper, err_est, err_code));
  return result;
  
  // ,
  // Rcpp::Named("error") = err_est
}

// // [[Rcpp::export]]
// double intensity_cdf_func(double upper, double mu,double beta_v11,double beta_v12,double alpha,double value ) {
// 
//   return (CDF_intensity(mu, beta_v11,beta_v12,alpha,0,upper)-value);
// 
// }


// [[Rcpp::export]]
double intensity_cdf_func(double upper, double mu,double beta_v11,double beta_v12,double alpha,double value ) {
  double nu=std::exp(beta_v11) ;
  double kappa=std::exp(beta_v12)+1 ;
  double cdf=R::pgamma(upper,kappa,nu/(kappa-1),1,0);
  double int_intensity=(std::exp(mu)*(upper)+alpha*cdf);
  double cdf_out= 1-std::exp(-int_intensity);
  return (cdf_out-value);
}



// [[Rcpp::export]]
double CDF_hazard(double shape, double beta_s, double e_y, double alpha, double beta_alpha, double h0, double beta_sd,
                  double beta_sd_cum, double di, double tox, double lambda_tox, double t_r, arma::vec beta_l, double b_il3, double lower,  double upper,
                  int Ncov_l) {
  
  haz_func function(shape,beta_s,e_y,alpha,beta_alpha,h0,beta_sd,beta_sd_cum,di,tox,lambda_tox,t_r, beta_l,b_il3,Ncov_l);
  double err_est;
  int err_code;
  const double result = 1-std::exp(-Numer::integrate(function, lower, upper, err_est, err_code));
  return result;
  
  // Rcpp::Named("error") = err_est
}

// [[Rcpp::export]]
double hazard_cdf_func(double t_r, double upper, double shape, double beta_s, double e_y, double alpha, double beta_alpha, double h0, 
                       double beta_sd, double beta_sd_cum, double di, double tox, double lambda_tox,arma::vec beta_l, double b_il3, double value, int Ncov_l ) {
  
  return (CDF_hazard(shape, beta_s,e_y,alpha,beta_alpha,h0, beta_sd,
                     beta_sd_cum, di, tox,lambda_tox,t_r,beta_l,b_il3,t_r,upper,Ncov_l)-value);
  
  // Rcpp::Named("error") = err_est
}


// [[Rcpp::export]]
double cum_haz_func(double t_r, double upper, double shape, double beta_s, double e_y, double alpha, double beta_alpha, double h0, 
                    double beta_sd, double beta_sd_cum, double di, double tox, double lambda_tox,arma::vec beta_l, double b_il3, double value, int Ncov_l ) {
  
  haz_func function(shape,beta_s,e_y,alpha,beta_alpha,h0,beta_sd,beta_sd_cum,di,tox,lambda_tox,t_r, beta_l,b_il3,Ncov_l);
  double err_est;
  int err_code;
  const double result = Numer::integrate(function, t_r, upper, err_est, err_code)-value;
  return result;
}

// [[Rcpp::export]]
double brents_fun_intensity(double lower, double upper, double tol, unsigned int max_iter, double mu,double beta_v11,double beta_v12,double alpha, double value) {
  double a = lower;
  double b = upper;
  double fa = intensity_cdf_func(a,mu,beta_v11,beta_v12,alpha,value);	// calculated now to save function calls
  double fb = intensity_cdf_func(b,mu,beta_v11,beta_v12,alpha,value);	// calculated now to save function calls
  double fs = 0;		// initialize 
  
  if (!(fa * fb < 0))
  {
    std::cout << "Intensity Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
    return -50;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < max_iter; ++iter)
  {
    // Rcout << iter <<std::endl;
    
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tol)
    {
      // std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
      return s;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tol) ) ||
         ( !mflag && (std::abs(c-d) < tol))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = intensity_cdf_func(s,mu,beta_v11,beta_v12,alpha,value);	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s) 
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for 
  
  std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
  return(-50);
} 

// [[Rcpp::export]]
double brents_fun_hazard(double lower, double upper, double tol, unsigned int max_iter, double shape, double beta_s, double e_y, double alpha,
                         double beta_alpha, double h0, double t_r,  double beta_sd,
                         double beta_sd_cum, double di, double tox, double lambda_tox, arma::vec beta_l, double b_il3, double value, int Ncov_l) {
  double a = lower;
  double b = upper;
  double fa = hazard_cdf_func(t_r, a, shape, beta_s, e_y, alpha, beta_alpha, h0, beta_sd,
                              beta_sd_cum, di, tox, lambda_tox, beta_l,  b_il3, value,Ncov_l);	// calculated now to save function calls
  double fb = hazard_cdf_func(t_r, b, shape, beta_s, e_y, alpha, beta_alpha, h0, beta_sd,
                              beta_sd_cum, di, tox,lambda_tox, beta_l, b_il3, value,Ncov_l);	// calculated now to save function calls
  double fs = 0;		// initialize 
  
  if (!(fa * fb < 0))
  {
    std::cout << "Hazard Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
    // std::cout << fa << std::endl; // throws exception if root isn't bracketed
    // std::cout << fb << std::endl; // throws exception if root isn't bracketed
    
    return -50;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < max_iter; ++iter)
  {
    // Rcout << iter <<std::endl;
    
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tol)
    {
      // std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
      return s;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tol) ) ||
         ( !mflag && (std::abs(c-d) < tol))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = hazard_cdf_func(t_r, s, shape, beta_s, e_y, alpha, beta_alpha, h0, beta_sd,
                         beta_sd_cum, di, tox, lambda_tox, beta_l, b_il3,value,Ncov_l );	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s) 
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for 
  
  std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
  return(-50);
    
} 


// [[Rcpp::export]]
double brents_fun_cum_haz(double lower, double upper, double tol, unsigned int max_iter, double shape, double beta_s, double e_y, double alpha,
                          double beta_alpha, double h0, double t_r,  double beta_sd,
                          double beta_sd_cum, double di, double tox, double lambda_tox, arma::vec beta_l, double b_il3, double value, int Ncov_l) {
  double a = lower;
  double b = upper;
  double fa = cum_haz_func(t_r, a, shape, beta_s, e_y, alpha, beta_alpha, h0, beta_sd,
                           beta_sd_cum, di, tox, lambda_tox, beta_l,b_il3, value,Ncov_l );	// calculated now to save function calls
  double fb = cum_haz_func(t_r, b, shape, beta_s, e_y, alpha, beta_alpha, h0, beta_sd,
                           beta_sd_cum, di, tox,lambda_tox, beta_l,b_il3,  value,Ncov_l );	// calculated now to save function calls
  double fs = 0;		// initialize 
  
  if (!(fa * fb < 0))
  {
    std::cout << "Hazard Signs of f(lower_bound) and f(upper_bound) must be opposites" << std::endl; // throws exception if root isn't bracketed
    // std::cout << fa << std::endl; // throws exception if root isn't bracketed
    // std::cout << fb << std::endl; // throws exception if root isn't bracketed
    
    return -50;
  }
  
  if (std::abs(fa) < std::abs(b))	// if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
  {
    std::swap(a,b);
    std::swap(fa,fb);
  }
  
  double c = a;			// c now equals the largest magnitude of the lower and upper bounds
  double fc = fa;			// precompute function evalutation for point c by assigning it the same value as fa
  bool mflag = true;		// boolean flag used to evaluate if statement later on
  double s = 0;			// Our Root that will be returned
  double d = 0;			// Only used if mflag is unset (mflag == false)
  
  for (unsigned int iter = 1; iter < max_iter; ++iter)
  {
    // Rcout << iter <<std::endl;
    
    // stop if converged on root or error is less than tolerance
    if (std::abs(b-a) < tol)
    {
      // std::cout << "After " << iter << " iterations the root is: " << s << std::endl;
      return s;
    } // end if
    
    if (fa != fc && fb != fc)
    {
      // use inverse quadratic interopolation
      s =	  ( a * fb * fc / ((fa - fb) * (fa - fc)) )
      + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
      + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }
    else
    {
      // secant method
      s = b - fb * (b - a) / (fb - fa);
    }
    
    // checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
    if (	( (s < (3 * a + b) * 0.25) || (s > b) ) ||
         ( mflag && (std::abs(s-b) >= (std::abs(b-c) * 0.5)) ) ||
         ( !mflag && (std::abs(s-b) >= (std::abs(c-d) * 0.5)) ) ||
         ( mflag && (std::abs(b-c) < tol) ) ||
         ( !mflag && (std::abs(c-d) < tol))	)
    {
      // bisection method
      s = (a+b)*0.5;
      
      mflag = true;
    }
    else
    {
      mflag = false;
    }
    
    fs = cum_haz_func(t_r, s, shape, beta_s, e_y, alpha, beta_alpha, h0, beta_sd,
                      beta_sd_cum, di, tox, lambda_tox, beta_l, b_il3, value,Ncov_l );	// calculate fs
    d = c;		// first time d is being used (wasnt used on first iteration because mflag was set)
    c = b;		// set c equal to upper bound
    fc = fb;	// set f(c) = f(b)
    
    if ( fa * fs < 0)	// fa and fs have opposite signs
    {
      b = s;
      fb = fs;	// set f(b) = f(s)
    }
    else
    {
      a = s;
      fa = fs;	// set f(a) = f(s) 
    }
    
    if (std::abs(fa) < std::abs(fb)) // if magnitude of fa is less than magnitude of fb
    {
      std::swap(a,b);		// swap a and b
      std::swap(fa,fb);	// make sure f(a) and f(b) are correct after swap
    }
    
  } // end for 
  
  std::cout<< "The solution does not converge or iterations are not sufficient" << std::endl;
  return(-50);
    
}



// [[Rcpp::export]]
double random_visit_times_Rcpp_test( double mu,double beta_v11,double beta_v12,double y_r, double t_r, arma::vec theta_a, double k) {
  double alpha=k/(1+exp(theta_a[0]+theta_a[1]*y_r));
  unsigned int max_iter=pow(10,4);
  double tol=pow(10,-8);
  double value=as<double>(Rcpp::runif(1));
  double lower=0;
  double upper=pow(10,9);
  // double gamma_k=std::tgamma(std::exp(beta_v12)+1) ;
  double vis_time=brents_fun_intensity(lower,upper, tol, max_iter, mu,beta_v11,beta_v12,alpha, value)+t_r;
  return(vis_time);
}

// [[Rcpp::export]]
double random_visit_times_Rcpp_test_fixed( double t_r,  double n_days) {

  double vis_time=n_days+t_r;
  return(vis_time);
}

// [[Rcpp::export]]
double toxicity_calc_RCpp(double d_i,double t_i,double t_r, double lambda_tox,double tox_recent) {
  // double tox=100*tox_recent*exp(-(t_i-t_r)/lambda_tox) + (1-exp(-(t_i-t_r)/lambda_tox ) )*lambda_tox*d_i;
  double tox=tox_recent*exp(-(t_i-t_r)/lambda_tox) + (1-exp(-(t_i-t_r)/lambda_tox ) )*d_i;
  
  return(std::max(tox,0.0));
}


// [[Rcpp::export]]
List random_hazard_times_Rcpp_test_exp( double t_max, double shape,double beta_s,double e_y, double beta_alpha,double h0, arma::vec theta_a,
                                        double y_r, double t_r, double beta_sd,double beta_sd_cum, 
                                        double di, double tox, double lambda_tox, double cum_haz_prev,arma::vec beta_l, double b_il3, int Ncov_l, double k ) {
  bool death;
  double alpha=k/(1+exp(theta_a[0]+theta_a[1]*y_r));
  unsigned int max_iter=pow(10,4);
  double tol=pow(10,-8);
  double value=std::log(2);
  double death_time=pow(10,20);
  double cum_haz=cum_haz_prev+integrate_hazard(shape,beta_s,e_y,alpha,beta_alpha,h0,
                                                          beta_sd,beta_sd_cum, di, tox,lambda_tox,t_r,beta_l,b_il3,
                                                          t_r, t_max,Ncov_l);  
  
  if (cum_haz>=value){
    death=true;
    // cum_haz=haz*(death_time-t_r)+cum_haz_prev;
    death_time=brents_fun_cum_haz( t_r, t_max, tol, max_iter, shape,  beta_s,  e_y,  alpha,  beta_alpha,  h0, t_r, beta_sd,
                                   beta_sd_cum, di, tox,lambda_tox,beta_l,b_il3, value-cum_haz_prev,Ncov_l);
    cum_haz=std::log(2);
  }else{
    death=false;
    
    
  }
  List output;
  output["cum_haz"]=cum_haz;
  output["death_time"]=death_time;
  output["death"]=death;
  
  return(output);
}


// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::vec mvrnormArma_single(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols; 
  arma::vec Y = arma::randn(ncols);
  return mu + arma::chol(sigma)*Y;
}

class int_func_grad_betav2_Rcpp: public Numer::Func {
private:
  double y_r;
  double beta_v11;
  double beta_v12;
  arma::vec theta_a;
  double k;
public:
  int_func_grad_betav2_Rcpp(double y_r_, double beta_v11_, 
           double beta_v12_, arma::vec theta_a_, double k_): y_r(y_r_), beta_v11(beta_v11_), 
           beta_v12(beta_v12_), theta_a(theta_a_), k(k_){}
  
  double operator()(const double& x) const {
    double nu=std::exp(beta_v11) ;
    double kappa=std::exp(beta_v12)+1 ;
    double exp_v2=std::exp(beta_v12);
    double alpha=k/(1+exp(theta_a[0]+theta_a[1]*y_r));
    double den_calc=R::dgamma(x,kappa,nu/(kappa-1),0);
    double tmp1=( (exp_v2*(beta_v12-std::log(nu))+kappa)*pow(x,exp_v2)*std::exp(-exp_v2*(x)/nu)+
      std::exp(beta_v12-exp_v2*(x)/nu)*std::log(x)*pow(x,exp_v2)-
      std::exp(beta_v12-exp_v2*(x)/nu)*pow(x,kappa)/nu );
    double tmp=0;
    if (den_calc>(pow(10,-50) )){

       tmp=tmp1*
        den_calc/(pow(x,exp_v2)*std::exp(-exp_v2*(x)/nu)) -
        den_calc*R::digamma(kappa);
    }
    double out=tmp*alpha;

    return (out);
    
  } 
}; 

// [[Rcpp::export]]
double integrate_grad_betav2_Rcpp(const double &y_r,const double &beta_v11,const double &beta_v12,arma::vec theta_a,double k,
                                        const double &lower, const double &upper) {
  int_func_grad_betav2_Rcpp function(y_r,beta_v11,beta_v12,theta_a,k);
  double err_est;
  int err_code;
  int subdiv=50;
  double eps_abs= pow(10,-2);
  double eps_rel= pow(10,-2);
  // double result = function(upper);
  double result = Numer::integrate(function, lower, upper, err_est, err_code, subdiv,
                                         eps_abs,eps_rel, Integrator<double>::GaussKronrod15);
  return result;
}

// [[Rcpp::export]]
double function_grad_betav2_Rcpp(const double &y_r,const double &beta_v11,const double &beta_v12,arma::vec theta_a,double k,
                                  const double &val) {
  int_func_grad_betav2_Rcpp function(y_r,beta_v11,beta_v12,theta_a,k);
  double result = function(val);
  return result;
}

// [[Rcpp::export]]
double beta_v2_grad_Rcpp(arma::vec J_cumsum,arma::vec J_cumsum_t, arma::vec J, 
                         int Npat, arma::vec Y, arma::vec Ts,double beta_v11, 
                    double beta_v12, double mu, arma::mat theta_as, 
                    arma::vec ks, arma::vec rewards) {
  double int_val=0;
  double sum_val=0;
  int min_it=0;
  int min_it_t=0;
  arma::vec theta_a(2);
  double k;
  double mean_reward=mean(rewards);
  for (int i=0; i<Npat ; ++i){
    theta_a=theta_as.col(i);
    k=ks(i);
    int max_it=J_cumsum[i]-1;
    int max_it_t=J_cumsum_t[i]-1;
    int j_i=J(i);
    arma::vec y_i=Y.subvec(min_it,max_it);
    arma::vec time=Ts.subvec(min_it_t,max_it_t);
    for (int j=1; j<(j_i+1) ; ++j){
      double x=time[j]-time[j-1];
      double y_r=y_i[j-1];
      double alpha=k/(1+exp(theta_a[0]+theta_a[1]*y_r));
      double int_ij_grad=function_grad_betav2_Rcpp(y_r, beta_v11,beta_v12, theta_a, k, x);
      double int_ij=test_intensity(x,mu,beta_v11, beta_v12, alpha);
      int_val=int_val-integrate_grad_betav2_Rcpp(y_r, beta_v11,beta_v12, theta_a, k, 0, x)*
        (rewards[i] - mean_reward) ;
      sum_val=sum_val+int_ij_grad/int_ij*(rewards[i] - mean_reward) ;
      
    }
    min_it=max_it+1;
    min_it_t=max_it_t+1;
  }
  double output= int_val+sum_val;
  return(output);
}


class int_func_grad_betav1_Rcpp: public Numer::Func {
private:
  double y_r;
  double beta_v11;
  double beta_v12;
  arma::vec theta_a;
  double k;
public:
  int_func_grad_betav1_Rcpp(double y_r_, double beta_v11_, 
                            double beta_v12_, arma::vec theta_a_, double k_): y_r(y_r_), beta_v11(beta_v11_), 
                            beta_v12(beta_v12_), theta_a(theta_a_), k(k_){}
  
  double operator()(const double& x) const {
    double nu=std::exp(beta_v11) ;
    double kappa=std::exp(beta_v12)+1 ;
    double alpha=k/(1+exp(theta_a[0]+theta_a[1]*y_r));
    double den_calc=R::dgamma(x,kappa,nu/(kappa-1),0);
    double intens=alpha*den_calc;
    double out=( (kappa-1)*(x)/nu-kappa )*intens;
    return (out);
    
  } 
}; 

// [[Rcpp::export]]
double integrate_grad_betav1_Rcpp(const double &y_r,const double &beta_v11,const double &beta_v12,arma::vec theta_a,double k,
                                  const double &lower, const double &upper) {
  int_func_grad_betav1_Rcpp function(y_r,beta_v11,beta_v12,theta_a,k);
  double err_est;
  int err_code;
  int subdiv=50;
  double eps_abs= pow(10,-2);
  double eps_rel= pow(10,-2);
  // double result = function(upper);
  double result = Numer::integrate(function, lower, upper, err_est, err_code, subdiv,
                                   eps_abs,eps_rel, Integrator<double>::GaussKronrod15);
  return result;
}

// [[Rcpp::export]]
double function_grad_betav1_Rcpp(const double &y_r,const double &beta_v11,const double &beta_v12,arma::vec theta_a,double k,
                                 const double &val) {
  int_func_grad_betav1_Rcpp function(y_r,beta_v11,beta_v12,theta_a,k);
  double result = function(val);
  return result;
}

// [[Rcpp::export]]
double beta_v1_grad_Rcpp(arma::vec J_cumsum,arma::vec J_cumsum_t, arma::vec J, 
                         int Npat, arma::vec Y, arma::vec Ts,double beta_v11, 
                         double beta_v12, double mu, arma::mat theta_as, 
                         arma::vec ks, arma::vec rewards) {
  double int_val=0;
  double sum_val=0;
  int min_it=0;
  int min_it_t=0;
  arma::vec theta_a(2);
  double k;
  double mean_reward=mean(rewards);
  for (int i=0; i<Npat ; ++i){
    theta_a=theta_as.col(i);
    k=ks(i);
    int max_it=J_cumsum[i]-1;
    int max_it_t=J_cumsum_t[i]-1;
    int j_i=J(i);
    arma::vec y_i=Y.subvec(min_it,max_it);
    arma::vec time=Ts.subvec(min_it_t,max_it_t);
    for (int j=1; j<(j_i+1) ; ++j){
      double x=time[j]-time[j-1];
      double y_r=y_i[j-1];
      double alpha=k/(1+exp(theta_a[0]+theta_a[1]*y_r));
      double int_ij_grad=function_grad_betav1_Rcpp(y_r, beta_v11,beta_v12, theta_a, k, x);
      double int_ij=test_intensity(x,mu,beta_v11, beta_v12, alpha);
      int_val=int_val-integrate_grad_betav1_Rcpp(y_r, beta_v11,beta_v12, theta_a, k, 0, x)*
        (rewards[i] - mean_reward) ;
      sum_val=sum_val+int_ij_grad/int_ij*(rewards[i] - mean_reward) ;
      
    }
    min_it=max_it+1;
    min_it_t=max_it_t+1;
  }
  double output= int_val+sum_val;
  return(output);
}




// [[Rcpp::export]]
Rcpp::List Data_Simu_Rcpp_single_pat_age_dgf( arma::mat beta_ls, arma::cube Bs, arma::vec sigma2_ls, 
                                              arma::vec beta_ss, arma::vec beta_sds,arma::vec beta_sd_cums,arma::vec lambda_toxs,
                                              arma::vec h0s,arma::vec beta_alphas,arma::vec shapes,
                                              int Npat,arma::vec beta_d, arma::vec sigma2_ds,double mu,double beta_v11, 
                                              double beta_v12,arma::mat theta_as,
                                              double beta_reward, arma::mat b_ils, arma::vec ks, arma::vec cov_in,
                                              double mean_y_init, int Ncov_l, int Ncov_d,arma::vec dosage_bounds){
  arma::vec beta_l (Ncov_l);
  arma::mat B(3,3);
  // arma::vec b_il(3);
  arma::vec cov_l(Ncov_l);
  arma::vec cov_d(Ncov_d);
  arma::vec cov_il(3);
  arma::vec theta_a(2);
  double sigma2_l;
  double k;
  double beta_s;
  double beta_sd;
  double beta_sd_cum;
  double h0;
  double beta_alpha;
  double shape;
  double tox;
  double death_i;
  double lambda_tox;
  double reward;
  double sigma2_d;
  // double mean_y_init=5;
  // double w_bmi=0.1;
  List death;
  int size=pow(10,7);
  arma::vec death_all(Npat);
  arma::vec reward_all(Npat);
  
  arma::vec t_i(size);
  arma::vec y_i(size);
  // arma::vec y_last_0(size);
  arma::vec d_i(size);
  arma::vec ey_i(size);
  arma::vec Toxicity(size);
  arma::vec mean_rand_i(size);
  arma::vec mean_fixed_i(size);
  arma::vec id(size);
  arma::vec id_t(size);
  arma::vec cum_haz(size);  ////////////////////
  
  arma::vec b_il_mu{0,0,0};
  arma::vec b_il(3);
  int tot_censor=0;
  int current_ind=0;
  int current_ind_t=0;
  bool out_bound=false;
  double first_vis=0;
  
  for (int i=0; i<Npat; ++i){
    beta_l=beta_ls.col(i);
    theta_a=theta_as.col(i);
    B= Bs.slice(i);
    sigma2_l=sigma2_ls(i);
    sigma2_d=sigma2_ds(i);
    beta_s=beta_ss(i);
    beta_sd=beta_sds(i);
    beta_sd_cum=beta_sd_cums(i);
    h0=h0s(i);
    beta_alpha=beta_alphas(i);
    shape=shapes(i);
    lambda_tox=lambda_toxs(i);
    k=ks(i);
    t_i[current_ind_t]=0; 
    bool alive_ind=true;
    cov_d(0)=1;
    cov_d(1)=mean_y_init;
    cov_d(arma::span(2,(Ncov_d-1) ))=cov_in;
    d_i[current_ind]=as<double>(Rcpp::rnorm(1,dot(cov_d,beta_d),pow(sigma2_d,0.5)));
    b_il=b_ils.col(i);
    mean_rand_i[current_ind]=b_il[0];
    mean_fixed_i[current_ind]=mean_y_init;
    y_i[current_ind]=mean_y_init;
    // y_last_0[current_ind_t]=y_i[current_ind];
    
    Toxicity[current_ind]=0;
    id[current_ind]=i+1;
    id_t[current_ind_t]=i+1;
    ey_i[current_ind]=mean_y_init;
    cum_haz[current_ind_t]=0; ////////////////////
    cum_haz[current_ind_t+1]=0;  ////////////////////
    tox=0;
    int j=0;
    double censor_time=pow(10,10);
    int censor=1;
    while(alive_ind && censor==1){
      ++j;
      ++current_ind;
      ++current_ind_t;
      t_i[current_ind_t]=random_visit_times_Rcpp_test(mu,beta_v11,beta_v12,y_i[current_ind-1], t_i[current_ind_t-1], theta_a,k);
      if (j==1){
        first_vis=t_i[current_ind_t];
      }
      if(t_i[current_ind_t]>censor_time){
        censor=0;
        tot_censor=tot_censor+1;
        
      } else{
        death=random_hazard_times_Rcpp_test_exp(t_i[current_ind_t],shape, beta_s, ey_i[current_ind-1], beta_alpha, h0, theta_a, ////////////////////////////////////////
                                                y_i[current_ind-1], t_i[current_ind_t-1], beta_sd, beta_sd_cum, d_i[current_ind-1], tox,lambda_tox,  cum_haz[current_ind_t-1],beta_l, b_il[2],Ncov_l,k);
        if (j>1){
          cum_haz[current_ind_t]=death["cum_haz"];    //////////////////////////////
        }
        if(death["death"] && j>1){
          alive_ind=false;
          id_t[current_ind_t]=i+1;
          // y_last_0[current_ind_t]=0;
          ++current_ind_t;
          death_i=as<double>(death["death_time"])-first_vis;
          // death_i=as<double>(death["death_time"]);
          death_all[i]=death_i;
          reward=log(death_i-beta_reward*(j+1));
          reward_all[i]=(reward);
        } else{
          // Rcout << "long begin" <<std::endl;
          id[current_ind]=i+1;
          id_t[current_ind_t]=i+1; 
          cov_l(0)=1;
          cov_l(1)=d_i[current_ind-1];
          cov_l(arma::span(2,(Ncov_l-3) ))=cov_in;
          cov_l(Ncov_l-2)=t_i[current_ind_t];
          cov_l(Ncov_l-1)=pow(t_i[current_ind_t],2);
          
          // cov_l={1,d_i[current_ind-1],age_i,DGF_i,diab_i,typedonor_i,BMI_i,t_i[current_ind_t],pow(t_i[current_ind_t],2)};
          cov_il={1,d_i[current_ind-1],t_i[current_ind_t]};
          // Rcout << "cov_y done" <<std::endl;
        
          mean_rand_i[current_ind]=dot(cov_il,b_il);
          mean_fixed_i[current_ind]=dot(cov_l,beta_l);
          ey_i[current_ind]=mean_rand_i[current_ind]+mean_fixed_i[current_ind];
          y_i[current_ind]=as<double>(Rcpp::rnorm(1,ey_i[current_ind],pow(sigma2_l,0.5)));
          // y_last_0[current_ind_t]=y_i[current_ind];
          // Rcout << "Y done" <<std::endl;
          cov_d(1)=y_i[current_ind];
          d_i[current_ind]=as<double>(Rcpp::rnorm(1,dot(cov_d,beta_d),pow(sigma2_d,0.5)));
          if (d_i[current_ind]<dosage_bounds(0) || d_i[current_ind]>dosage_bounds(1)){
            
            out_bound=true;
          }
          tox=toxicity_calc_RCpp( d_i[current_ind-1], t_i[current_ind_t], t_i[current_ind_t-1],  lambda_tox, tox);
          // tox=1;
          Toxicity[current_ind]=tox;
          if (j==1){
            Toxicity[current_ind]=0;
            tox=0;
          }
          // Rcout << "Long done" <<std::endl;
          
        }
        
      }
      
    }
  }
  List output;
  output["Ts"]=t_i.head(current_ind_t);
  output["Y"]=y_i.head(current_ind);
  // output["Y_last_0"]=y_last_0.head(current_ind_t);
  output["D"]=d_i.head(current_ind);
  // output["BMI"]=BMI_i.head(current_ind);
  output["EY"]=ey_i.head(current_ind);
  // output["b_il"]=b_il;
  output["Toxicity"]=Toxicity.head(current_ind);
  output["id"]=id.head(current_ind);
  output["id_t"]=id_t.head(current_ind_t);
  // output["X0"]=Age.head(current_ind);
  // output["DGF"]=DGF.head(current_ind);
  // output["diab"]=diab.head(current_ind);
  // output["typedonor"]=typedonor.head(current_ind);
  // output["BMI"]=BMI.head(current_ind);
  output["Npat"]=Npat;
  output["death_time_i"]=death_all;
  output["reward_i"]=reward_all;
  output["cum_haz"]=cum_haz.head(current_ind_t);
  output["out_bound"]=out_bound;
  
  return(output);
  
}


// [[Rcpp::export]]
Rcpp::List Data_Simu_Rcpp_single_pat_fixed_visits( arma::mat beta_ls, arma::cube Bs, arma::vec sigma2_ls, 
                                              arma::vec beta_ss, arma::vec beta_sds,arma::vec beta_sd_cums,arma::vec lambda_toxs,
                                              arma::vec h0s,arma::vec beta_alphas,arma::vec shapes,
                                              int Npat,arma::vec beta_d, arma::vec sigma2_ds,double mu,double beta_v11, 
                                              double beta_v12,arma::mat theta_as,
                                              double beta_reward, arma::mat b_ils, arma::vec ks, arma::vec cov_in,
                                              double mean_y_init, int Ncov_l, int Ncov_d,arma::vec dosage_bounds, double n_days){
  arma::vec beta_l (Ncov_l);
  arma::mat B(3,3);
  arma::vec cov_l(Ncov_l);
  arma::vec cov_d(Ncov_d);
  arma::vec cov_il(3);
  arma::vec theta_a(2);
  double sigma2_l;
  double k;
  double beta_s;
  double beta_sd;
  double beta_sd_cum;
  double h0;
  double beta_alpha;
  double shape;
  double tox;
  double death_i;
  double lambda_tox;
  double reward;
  double sigma2_d;
  // double mean_y_init=5;
  // double w_bmi=0.1;
  List death;
  int size=pow(10,7);
  arma::vec death_all(Npat);
  arma::vec reward_all(Npat);
  
  arma::vec t_i(size);
  arma::vec y_i(size);
  
  // arma::vec y_last_0(size);
  
  arma::vec d_i(size);
  arma::vec ey_i(size);
  arma::vec Toxicity(size);
  arma::vec mean_rand_i(size);
  arma::vec mean_fixed_i(size);
  arma::vec id(size);
  arma::vec id_t(size);
  arma::vec cum_haz(size);  ////////////////////
  
  arma::vec b_il_mu{0,0,0};
  arma::vec b_il(3);
  int tot_censor=0;
  int current_ind=0;
  int current_ind_t=0;
  bool out_bound=false;
  double first_vis=0;
  
  for (int i=0; i<Npat; ++i){
    beta_l=beta_ls.col(i);
    theta_a=theta_as.col(i);
    B= Bs.slice(i);
    sigma2_l=sigma2_ls(i);
    sigma2_d=sigma2_ds(i);
    beta_s=beta_ss(i);
    beta_sd=beta_sds(i);
    beta_sd_cum=beta_sd_cums(i);
    h0=h0s(i);
    beta_alpha=beta_alphas(i);
    shape=shapes(i);
    lambda_tox=lambda_toxs(i);
    k=ks(i);
    t_i[current_ind_t]=0; 
    bool alive_ind=true;
    cov_d(0)=1;
    cov_d(1)=mean_y_init;
    cov_d(arma::span(2,(Ncov_d-1) ))=cov_in;
    d_i[current_ind]=as<double>(Rcpp::rnorm(1,dot(cov_d,beta_d),pow(sigma2_d,0.5)));
    b_il=b_ils.col(i);
    mean_rand_i[current_ind]=b_il[0];
    mean_fixed_i[current_ind]=mean_y_init;
    y_i[current_ind]=mean_y_init;
    // y_last_0[current_ind_t]=y_i[current_ind];
    
    Toxicity[current_ind]=0;
    id[current_ind]=i+1;
    id_t[current_ind_t]=i+1;
    ey_i[current_ind]=mean_y_init;
    cum_haz[current_ind_t]=0; ////////////////////
    cum_haz[current_ind_t+1]=0;  ////////////////////
    tox=0;
    int j=0;
    double censor_time=pow(10,10);
    int censor=1;
    while(alive_ind && censor==1){
      ++j;
      // Rcout << cum_haz[current_ind_t] <<std::endl;
      ++current_ind;
      ++current_ind_t;
      // t_i[current_ind_t]=random_visit_times_Rcpp_test(mu,beta_v11,beta_v12,y_i[current_ind-1], t_i[current_ind_t-1], theta_a,k);
      t_i[current_ind_t]=random_visit_times_Rcpp_test_fixed( t_i[current_ind_t-1],  n_days);
      if (j==1){
        first_vis=t_i[current_ind_t];
      }
      // Rcout << "Visit done" <<std::endl;
      if(t_i[current_ind_t]>censor_time){
        censor=0;
        tot_censor=tot_censor+1;
        
      } else{
        death=random_hazard_times_Rcpp_test_exp(t_i[current_ind_t],shape, beta_s, ey_i[current_ind-1], beta_alpha, h0, theta_a, ////////////////////////////////////////
                                                y_i[current_ind-1], t_i[current_ind_t-1], beta_sd, beta_sd_cum, d_i[current_ind-1], tox,lambda_tox,  cum_haz[current_ind_t-1],beta_l, b_il[2],Ncov_l,k);
        // Rcout << "death done" <<std::endl;
        if (j>1){
          cum_haz[current_ind_t]=death["cum_haz"];    //////////////////////////////
        }
        if(death["death"] && j>1){
          alive_ind=false;
          id_t[current_ind_t]=i+1;
          // y_last_0[current_ind_t]=0;
          ++current_ind_t;
          death_i=as<double>(death["death_time"])-first_vis;
          // death_i=as<double>(death["death_time"]);
          death_all[i]=death_i;
          reward=log(death_i-beta_reward*(j+1));
          reward_all[i]=(reward);
        } else{
          // Rcout << "long begin" <<std::endl;
          id[current_ind]=i+1;
          id_t[current_ind_t]=i+1; 
          cov_l(0)=1;
          cov_l(1)=d_i[current_ind-1];
          cov_l(arma::span(2,(Ncov_l-3) ))=cov_in;
          cov_l(Ncov_l-2)=t_i[current_ind_t];
          cov_l(Ncov_l-1)=pow(t_i[current_ind_t],2);
          
          // cov_l={1,d_i[current_ind-1],age_i,DGF_i,diab_i,typedonor_i,BMI_i,t_i[current_ind_t],pow(t_i[current_ind_t],2)};
          cov_il={1,d_i[current_ind-1],t_i[current_ind_t]};
          // Rcout << "cov_y done" <<std::endl;
          
          mean_rand_i[current_ind]=dot(cov_il,b_il);
          mean_fixed_i[current_ind]=dot(cov_l,beta_l);
          ey_i[current_ind]=mean_rand_i[current_ind]+mean_fixed_i[current_ind];
          y_i[current_ind]=as<double>(Rcpp::rnorm(1,ey_i[current_ind],pow(sigma2_l,0.5)));
          // y_last_0[current_ind_t]=y_i[current_ind];
          // Rcout << "Y done" <<std::endl;
          cov_d(1)=y_i[current_ind];
          d_i[current_ind]=as<double>(Rcpp::rnorm(1,dot(cov_d,beta_d),pow(sigma2_d,0.5)));
          // d_i[current_ind]=std::max(as<double>(Rcpp::rnorm(1,dot(cov_d,beta_d),pow(sigma2_d,0.5))),0.0 ) ;
          // if (d_i[current_ind]<0 || d_i[current_ind]>2.75){
          if (d_i[current_ind]<dosage_bounds(0) || d_i[current_ind]>dosage_bounds(1)){
            
            out_bound=true;
          }
          tox=toxicity_calc_RCpp( d_i[current_ind-1], t_i[current_ind_t], t_i[current_ind_t-1],  lambda_tox, tox);
          // tox=1;
          Toxicity[current_ind]=tox;
          if (j==1){
            Toxicity[current_ind]=0;
            tox=0;
          }
          // Rcout << "Long done" <<std::endl;
          
        }
        
      }
      
    }
  }
  List output;
  output["Ts"]=t_i.head(current_ind_t);
  output["Y"]=y_i.head(current_ind);
  // output["Y_last_0"]=y_last_0.head(current_ind_t);
  output["D"]=d_i.head(current_ind);
  output["EY"]=ey_i.head(current_ind);
  output["Toxicity"]=Toxicity.head(current_ind);
  output["id"]=id.head(current_ind);
  output["id_t"]=id_t.head(current_ind_t);
  output["Npat"]=Npat;
  output["death_time_i"]=death_all;
  output["reward_i"]=reward_all;
  output["cum_haz"]=cum_haz.head(current_ind_t);
  output["out_bound"]=out_bound;
  
  return(output);
  
}

