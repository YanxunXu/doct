# doct
"doct" stands for Decisions Optimized in Continuous Time. 

This R package implements the method from the paper "Personalized Dynamic Treatment Regimes in Continuous Time: A Bayesian Joint Model for Optimizing Clinical Decisions with Timing. "

When installing the R package using "doct_1.0.tar.gz", please use the following script to first install the dependent packages if they are not installed on your computer. 

dep_packages <- c("Rcpp","RcppArmadillo", "RcppEigen", "RcppNumerical", "MASS", "pracma", "mvtnorm", "LaplacesDemon","stats", "GoFKernel","survival")
new.packages <- dep_packages[!(dep_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

"simulation.Rdata" is a simulated dataset mimicking the real dataset. The detailed description on how to generate the simulated dataset can be found in the book chapter "Modeling and Optimizing Dynamic Treatment Regimens in Continuous Time. "







