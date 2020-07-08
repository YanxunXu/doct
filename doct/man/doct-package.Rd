\name{doct-package}
\alias{doct-package}
\alias{doct}
\docType{package}
\title{
\packageTitle{doct}
}
\description{
\packageDescription{doct}
}
\details{
The main functions in the package are mcmc_joint and SGD_run. mcmc_joint runs the Bayesian parameter estimation and collects the MCMC samples
to be used SGD_run, which optimizes the dosage and visitation scheduling policy. The other functions in this package are for analyzing the results of these two main functions. For example, DIC_calculation outputs the DIC for simulated data analysis. This package accompanies our statistical methodology paper.

Our joint model has three submodels. For patient i, the visitation intensity function and dosing submodel is expressed below for \eqn{t \in (t_{i,j} , t_{i, j+1}]} :
\deqn{  \lambda_{i}(t)=\exp(\mu) + \ \ \alpha_{i,j} (t-t_{i,j})^{\kappa-1} e^{-\gamma(t-t_{i,j})} \frac{\gamma^{\kappa}}{\Gamma(\kappa)}}
\deqn{\alpha_{i,j}= \frac{\xi}{1+\exp((1,y_{i,j}) \bm{\beta_\alpha})} \ \ \frac{\kappa-1}{\gamma}=\exp(\nu_1) , \ \ \kappa=\exp(\nu_2)+1 }
\deqn{ d_{i,j}=  (1,y_{i,j}, \bm{x}_{i} ) \bm{\beta}_{d}  + \epsilon_d, \ \ \epsilon_d \sim \mathrm{Normal}(0,\sigma_{d}^2)}

The longitudinal submodel for the log creatinine levels is represented below, with the following likelihood function:
\deqn{ y_i(t) = y_i^*(t) + \epsilon_l = \bm{z}_i(t) \bm{\beta}_{l}   + \bm{r}_i(t) \bm{b}_i+ \epsilon_l , \ \ \epsilon_l \sim \mathrm{Normal}(0,\sigma_{l}^2),\ \ \bm{b}_i \sim \mathrm{Normal}( \bf{0}, \bm{\Sigma}_b)}
\deqn{	\bm{z}_i(t)=( 1,d_i(t), \bm{x}_{i}, t,t^2) , \ \  \bm{r}_i(t)=(1,d_i(t),t), }

The survival submodel is represented using the hazard function:
\deqn{ h_{i}(t)=\exp\Big(-(\beta_{s1} y_i^*(t) + \beta_{s2} d_{i}(t) + \beta_{s3} \mathrm{Tox}_i(t)  + \beta_{s4} \alpha_{i}(t) +h_0 )\Big)st^{s-1} }

}
\author{
\packageAuthor{doct}

Maintainer: \packageMaintainer{doct}
}
