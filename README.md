# Approximate Bayesian Compuataion(ABC)
ABC is a powerful tools to estimates parameters in mathematical models.  
## Summary Statistic
Summary statistics is a crucial component in ABC. It is used to estimate the similarity between observated data and simulated data. 
As summary statistics is used to replace the likelihood, the selection of summary statistics has one important criteria, sensitivity.  
A sensitive summary statistics can reflect the difference between two dataset well enough.  

## ABC Rejection
ABC rejection is the easiest among all the algorithms in ABC. However, the efficiency is relative low.  
We first draw the sample from proposal distribution, then use the sample to simulate data. 
The similarity of two data set is calculated by the summary statistic.  These steps obtain a set of particles { $\theta_i$, $\rho_i$} for i in 1, ..., N  
Sorted the data by $\rho$ and drop the sample from largest $\rho$ until the it reach the threshold that given in adnvance.
ABC rejection is commonly used to select the summary statistics. 

## ABC MCMC
ABC MCMC use the concept of MCMC to filter the samples.  
The samples are still draw from proposal distribution and use to simulate the data. The similarity of simulated data and observed is still calculated by summary statistics to replace the likelihood.  
Apply the Metropolis-Hasting ratio to calculate the acceptance rate. 
