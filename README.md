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
Algorithm 1: ABC Rejection  
1	Define threshold, ϵ, target number of particle, T, initialise t = 1.  
2	While t < T do  
3		Draw μ_t,δ_t,ν_t from unif (0, 1)  
  
4		Simulate data, y, from Gillespie Algorithm   
5		Calculate the closeness of Simulated data and observation, ϵ_t= ρ(y,x)  
  
6		If ϵ_t < threshold:  
7			store the particle {μ_t,δ_t,ν_t,ϵ_t}  
  
8			t = t + 1  

## ABC MCMC
ABC MCMC use the concept of MCMC to filter the samples.  
The samples are still draw from proposal distribution and use to simulate the data. The similarity of simulated data and observed is still calculated by summary statistics to replace the likelihood.  
Apply the Metropolis-Hasting ratio to calculate the acceptance rate.  
Algorithm 2: ABC MCMC  
1	Draw μ_0,δ_0,ν_0 from unif(0, 1) and simulate x_0  
  
2	Calculus closeness ϵ_0=ρ(y,x_0)  
  
3	For t = 1 to N do  
4		Draw candidate μ^*,δ^*,ν^* from q(.|μ_(t-1),δ_(t-1),ν_(t-1))  
  
5		Simulate x_t according to μ^*,δ^*,ν^*  
  
6		Calculus closeness ψ_t=ρ(y,x_t)  
  
7		Compute MH ratio r=  (π(μ^* )π(δ^* )π(ν^* ) ψ^* q(μ_(t-1),δ_(t-1),ν_(t-1) |μ^*,δ^*,ν^*))/(π(μ_(t-1) )π(δ_(t-1) )π(ν_(t-1) ) ψ^(t-1)  q(μ^*,δ^*,ν^* |μ_(t-1),δ_(t-1),ν_(t-1)))  
  
8		If U(0,1) < r then  
  
9			Accept candidate and ψ^t= ψ^*  
10		else  
11			μ_t=μ_(t-1),δ_t= δ_(t-1),ν_t=ν_(t-1),ψ^t=ψ^(t-1)  

## AMC SMC
Algorithm 3: SMC ABC  
1	Set α, the percentage of particles kept each iteration. Define target closeness (tolerance), ϵ_T  
2	Draw N particles from prior, simulate data and calculus the closeness. This obtains {μ_t,δ_t,ν_t 〖,ρ_t}〗_(t=1)^N.  
  
3	While ϵ_N<ϵ_T do  
4		Sort the particles increasingly by ρ_t , set tolerance, ϵ_t= ρ_(N-αN)  
5		Select α*N particles with lowest closeness  
6		Compute the tuning parameters of the MCMC kernel q_t (.|.)  
7		For j=N-αN+1 to N do  
8			For k = 1 to S_t do  
9				Propose θ^* from q_t (.|θ^j) and simulate data x according to θ^*  
10				Compute MH ratio r=  (π(θ^*)q(θ^j |θ^*))/(π(θ^j)q(θ^* |θ^j))1(ρ(x,y)<ϵ_t)  
11				If unif((0,1) < r do  
12					θ^j=θ^* and ρ^j= ρ(x,y)  
  
13		If MCMC acceptance rate is too small   
14			break  
