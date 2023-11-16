library(smfsb)
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

# least_squared_error = function(observation, simulated){
#   coral = lm(observation$coral ~ simulated$coral)
#   starfish = lm(observation$starfish ~ simulated$starfish)
#   x = (sum((data$coral - coral$fitted.values)^2))^(1/2)
#   y = (sum((data$starfish - starfish$fitted.values)^2))^(1/2)
#   return(x+y)
# }

# theta = (mu, delta, nu)
N = 10^5
mu = 0.5
delta = 0.5
nu = 0.5
phi = Euclidean_distance(data, sim_Gillespie(mu, delta, nu))
theta = data.frame((matrix(ncol = 4, nrow = N)))
colnames(theta) = c('mu', 'delta', 'nu', 'phi')
theta[1,] = c(mu, delta, nu, phi)
sd = diag(0.1, 3, 3)


for (i in 2:N){
  # Propose candidate parameters1
  candidate = mvrnorm(1, c(theta$mu[i-1], theta$delta[i-1], theta$nu[i-1]), sd)
  if (sum(candidate < 0) > 0 | sum(candidate > 1) > 0){
    theta[i, ] = theta[i - 1, ]
    next
  }
  # Simulate new data
  new_mu = candidate[1]
  new_delta = candidate[2]
  new_nu = candidate[3]
  simulated = sim_Gillespie(new_mu, new_delta, new_nu, C0 = 34, S0= 16,
                            times = seq(from = 0,to = 20, by = 2))
  
  # The summary statistic needs to be in contrast to ABC rejection. 
  # Least square error does not available
  phi = Euclidean_distance(data, simulated)
  # How to calculate r? 
  r = phi/theta$phi[i-1]
  if (runif(1) < r){
    update = update + 1
    theta[i, ] = c(new_mu, new_delta, new_nu, phi)
  }else{
    theta[i, ] = theta[i - 1, ]
  }
  
}
# theta = theta[order(theta$phi), ]
plot(density(theta$mu),main ='mu', xlab = 'mu')
plot(density(theta$delta),main ='delta', xlab = 'delta')
plot(density(theta$nu),main ='nu', xlab = 'nu')

