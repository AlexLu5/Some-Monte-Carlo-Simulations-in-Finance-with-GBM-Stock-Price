S0 <- 1
mu <- 0.25
sigma <- 0.5
r <- 0.15
K <- 1

T <- 1
M <- 10000
N <- 10000
set.seed(0)

t <- rep(NA,N*T+1)
Z <- matrix(rnorm(M*N*T),M,N*T)
W <- matrix(NA,M,N*T+1)
S <- matrix(NA,M,N*T+1)
SQ <- matrix(NA,M,N*T+1)
C <- matrix(NA,M,N*T+1)
P <- matrix(NA,M,N*T+1)

E <- matrix(NA,2,N*T+1)
Var <- matrix(NA,2,N*T+1)
d <- matrix(NA,2,N*T)
Call <- matrix(NA,2,N*T+1)
Put <- matrix(NA,2,N*T+1)

W[ ,1] <- 0
S[ ,1] <- S0
SQ[ ,1] <- S0
Call[2,1] <- max(S0-K,0)
Put[2,1] <- max(K-S0,0)

for (n in 1:(N*T+1)){
  t[n] <- (n-1)/N
  
  for (m in 1:M){
    if (n>=2){
      W[m,n] <- W[m,n-1]+1/sqrt(N)*Z[m,n-1]
      S[m,n] <- S0*exp((mu-sigma^2/2)*t[n]+sigma*W[m,n])
      SQ[m,n] <- S0*exp((r-sigma^2/2)*t[n]+sigma*W[m,n])
    }
    C[m,n] <- exp(-r*t[n])*max(SQ[m,n]-K,0)
    P[m,n] <- exp(-r*t[n])*max(K-SQ[m,n],0)
  }
  E[1,n] <- mean(S[ ,n])
  E[2,n] <- S0*exp(mu*t[n])
  Var[1,n] <- var(S[ ,n])*(M-1)/M
  Var[2,n] <- S0^2*exp(2*mu*t[n])*(exp(sigma^2*t[n])-1)
  
  Call[1,n] <- mean(C[ ,n])
  Put[1,n] <- mean(P[ ,n])
  
  if (n>=2){
    d[1,n-1] <- (log(S0)-log(K))/(sigma*sqrt(t[n]))+r*sqrt(t[n])/sigma+sigma*sqrt(t[n])/2
    d[2,n-1] <- (log(S0)-log(K))/(sigma*sqrt(t[n]))+r*sqrt(t[n])/sigma-sigma*sqrt(t[n])/2
    Call[2,n] <- S0*pnorm(d[1,n-1])-K*exp(-r*t[n])*pnorm(d[2,n-1])
    Put[2,n] <- K*exp(-r*t[n])*pnorm(-d[2,n-1])-S0*pnorm(-d[1,n-1])
  }
}
png(filename="True and Estimated Mean of Stock Price.png")
plot(t,E[1, ],type="l",col="blue",xlab="t", ylab="Stock Price Means", main="True and Estimated Mean of Stock Price")
lines(t,E[2, ],col="red")
legend('topleft', legend = c('Estimated', 'True'), col = c('blue', 'red'), lty=1)
dev.off()

png(filename="True and Estimated Variance of Stock Price.png")
plot(t,Var[1, ],type="l",col="blue",xlab="t", ylab="Stock Price Variances", main="True and Estimated Variance of Stock Price")
lines(t,Var[2, ],col="red")
legend('topleft', legend = c('Estimated', 'True'), col = c('blue', 'red'), lty=1)
dev.off()

png(filename="True and Estimated Call Price.png")
plot(t,Call[1, ],type="l",col="blue",xlab="t", ylab="Call Prices", main="True and Estimated Call Price")
lines(t,Call[2, ],col="red")
legend('topleft', legend = c('Estimated', 'True'), col = c('blue', 'red'), lty=1)
dev.off()

png(filename="True and Estimated Put Price.png")
plot(t,Put[1, ],type="l",col="blue",xlab="t", ylab="Put Prices", main="True and Estimated Put Price")
lines(t,Put[2, ],col="red")
legend('topleft', legend = c('Estimated', 'True'), col = c('blue', 'red'), lty=1)
dev.off()

print(paste0("Plots in ",getwd()))