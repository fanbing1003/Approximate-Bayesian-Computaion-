#install.packages('smfsb')
library(smfsb)
library(coda)
library(MASS)
data = read.csv('data_project.csv')
set.seed(13)

# mu: Coral + 1
# delta: Coral - 1, Starfish + 1
# nu: Starfish -1 
sim_Gillespie <- function(mu, delta, nu, C0 = 34, S0= 16,
                          times = seq(from = 0,to = 20, by = 2)){
  x = c(C0, S0)
  states = matrix(0, nrow = length(times), length(x))
  states[1,] <- x
  
  for( i in 2:length(times)){
    x <- smfsb::stepLVc(x0 = x, t0 = times[i-1], deltat = times[i] - times[i-1], th = c(mu, delta, nu))
    states[i, ] <- x
  }
  return(data.frame(times=times, coral = states[,1], starfish = states[,2]))
}

Euclidean_distance = function(observation, simulated){
  d = sum((observation$coral - simulated$coral)^2) + sum((observation$starfish - simulated$starfish)^2)
  return(d^(1/2))
}

# Drop half of the particle
alpha = 0.5
# Number of iteration on each particle to use for each move step
st = 100
# Draw N particle from prior, simulate data and calculate the Euclidean distance.
N = 10000
theta = data.frame((matrix(ncol = 4, nrow = N)))
colnames(theta) = c('mu', 'delta', 'nu', 'phi')
theta$mu = runif(N)
theta$delta = runif(N)
theta$nu = runif(N)
for (i in 1:N){
  simulated =  sim_Gillespie(theta$mu[i], theta$delta[i], theta$nu[i])
  theta$phi[i] = Euclidean_distance(data, simulated)
}


# Sort the particle set by phi
theta = theta[order(theta$phi), ]

# Target Tolerance
epslion_Target = 170

epslion_max = theta$phi[N]
accep = 0
resample_num = 5001

while (epslion_max > epslion_Target){
  #Sort the particle by phi, set up epsilon_max
  theta = theta[order(theta$phi), ]
  epslion_t = theta$phi[N - alpha*N]
  
  # set up for MCMC kernel 
  mean_mu = mean(theta$mu[N - alpha*N])
  mean_delta = mean(theta$delta[N - alpha*N])
  mean_nu = mean(theta$nu[N - alpha*N])
  sd_mu = sd(theta$mu[1:(alpha*N)])
  sd_delta = sd(theta$delta[1:(alpha*N)])
  sd_nu = sd(theta$nu[1:(alpha*N)])
  cov = matrix(c(sd_mu, 0, 0, 
                 0, sd_nu, 0, 
                 0, 0, sd_nu), 
               nrow = 3, ncol = 3, byrow = TRUE)
  

  for (j in resample_num:N){
    for (k in 1:1000){
      #  Simple RW for re-sample
      params = mvrnorm(1, c(mean_mu, mean_delta, mean_nu), cov)
      if (sum(params < 0) > 0 | sum(params > 1) > 0){
        next
      }
      simulated = sim_Gillespie(params[1], params[2], params[3])
      # Calculated MH acceptance rate
      if (Euclidean_distance(data, simulated) < epslion_t){
        theta$mu[j] = params[1]                
        theta$delta[j] = params[2]
        theta$nu[j] = params[3]
        theta$phi[j] = Euclidean_distance(data, simulated)
        accep = accep + 1
        break
      }
    }
  }
  #Rt = log(0.01)/log(1-accp/5000)
  par(mfrow = c(1, 3))
  print(sprintf('Current Tolerance: %f', epslion_t))
  print(sprintf('MCMC Accepted New Sample Number %f', accep))
  plot(density(theta$mu),main ='mu')
  plot(density(theta$delta),main ='delta')
  plot(density(theta$nu),main ='nu')
  theta = theta[order(theta$phi), ]
  epslion_max = theta$phi[N]
  if (epslion_max < epslion_Target){
    break
  }
  if (accep < 1000){
    break
  }
  accep = 0
}


