#install.packages('smfsb')
library(smfsb)
data = read.csv('data_project.csv')

# Mathematical model
# C(t) is the number of coral alive at time t
# S(t) is the number of crown of thorn starfish alive at time t

# This is the function that used to simulate data based on given parameter theta = (mu, delta, nu)
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
# Euclidean Distance
Euclidean_distance = function(observation, simulated){
  d = sum((observation$coral - simulated$coral)^2) + sum((observation$starfish - simulated$starfish)^2)
  return(d^(1/2))
}

least_squared_error = function(observation, simulated){
  coral = lm(observation$coral ~ simulated$coral)
  starfish = lm(observation$starfish ~ simulated$starfish)
  x = (sum((data$coral - coral$fitted.values)^2))^(1/2)
  y = (sum((data$starfish - starfish$fitted.values)^2))^(1/2)
  return(x+y)
}

# ABC Rejection
N = 10^5
mu = 0.5
delta = 0.5
nu = 0.5
theta = data.frame((matrix(ncol = 4, nrow = N)))
colnames(theta) = c('mu', 'delta', 'nu', 'phi')
theta[1,] = c(mu, delta, nu, 10000)

for (i in 2:N){
    new_mu = runif(1)
    new_delta = runif(1)
    new_nu = runif(1)
    simulated = sim_Gillespie(new_mu, new_delta, new_nu, C0 = 34, S0= 16,
                                times = seq(from = 0,to = 20, by = 2))
    theta[i, ] = c(new_mu, new_delta, new_nu, Euclidean_distance(data, simulated))
}

alpha = 0.01
threshold = quantile(theta$phi, alpha, na.rm=TRUE)
print(threshold)
row = 1
theta_final = data.frame((matrix(ncol = 4, nrow = N*alpha)))
colnames(theta_final) = c('mu', 'delta', 'nu', 'phi')

for (i in 1:N){
  if (theta$phi[i] <= threshold[[1]]){
    theta_final[row,] = theta[i, ]
    row = row + 1
  }
}
theta = theta[order(theta$phi), ]
theta$phi[1000]
theta$phi[500]
theta$phi[100]
# Mu
plot(density(theta$mu[1:100]), main = 'mu', col = 'blue', xlim = c(0,1), lty = 1, xlab = 'mu', lwd = 2)
lines(density(theta$mu[1:500]), main = 'mu', col = 'black', xlim = c(0,1), lty = 2, lwd = 2)
lines(density(theta$mu[1:1000]), main = 'mu', col = 'red', xlim = c(0,1), lty = 3, lwd = 2)
legend("topleft", legend = c('1% (224.5)', '0.5% (208.9)', '0.1% (158.1)'), col = c('red', 'black', 'blue'), lty=1:3, cex=0.8)


# Delta
plot(density(theta$delta[1:100]), main = 'delta', col = 'blue', xlim = c(0,0.2), lty = 1, xlab = 'delta', lwd = 2)
lines(density(theta$delta[1:500]), main = 'delta', col = 'black', xlim = c(0,0.2), lty = 2, lwd = 2)
lines(density(theta$delta[1:1000]), main = 'delta', col = 'red', xlim = c(0,0.2), lty = 3, lwd = 2)
legend("topright", legend = c('1% (224.5)', '0.5% (208.9)', '0.1% (158.1)'), col = c('red', 'black', 'blue'), lty=1:3, cex=0.8)

# nu
plot(density(theta$nu[1:100]), main = 'nu', col = 'blue', xlim = c(0,1), lty = 2, lwd = 2)
lines(density(theta$nu[1:1000]), main = 'nu', col = 'red', xlim = c(0,1), lty = 1, xlab = 'nu', lwd = 2)
lines(density(theta$nu[1:500]), main = 'nu', col = 'black', xlim = c(0,1), lty = 3, lwd = 2)
legend("topleft", legend = c('1% (224.5)', '0.5% (208.9)', '0.1% (158.1)'), col = c('red', 'black', 'blue'), lty=1:3, cex=0.8)

plot(density(theta$phi), xlim = c(50, 10000), main = 'Density of Rho (Euclidean Distance)', xlab = 'Rho')



