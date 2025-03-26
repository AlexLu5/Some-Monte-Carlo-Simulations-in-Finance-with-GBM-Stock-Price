using Distributions
using Plots
using Random
using Statistics

S0 = 1
mu = 0.25
sigma = 0.5
r = 0.15
K = 1

T = 1
N = 10000
M = 10000
Random.seed!(0)

t = [(n-1)/N for n in 1:N*T+1]
Z = randn(M,N*T+1)
W = fill(0.0,(M,N*T+1))
S = fill(0.0,(M,N*T+1))
SQ = fill(0.0,(M,N*T+1))
C = fill(0.0,(M,N*T+1))
P = fill(0.0,(M,N*T+1))

E = fill(0.0,(2,N*T+1))
Var = fill(0.0,(2,N*T+1))
d = fill(0.0,(2,N*T+1))
Call = fill(0.0,(2,N*T+1))
Put = fill(0.0,(2,N*T+1))

S[:,1] .= S0
SQ[:,1] .= S0
Call[2,1] = max(S0-K,0)
Put[2,1] = max(K-S0,0)

for n in 1:N*T+1
    for m in 1:M
        if n>=2
            W[m,n] = W[m,n-1]+1/sqrt(N)*Z[m,n-1]
            S[m,n] = S0*exp((mu-sigma^2/2)*t[n]+sigma*W[m,n])
            SQ[m,n] = S0*exp((r-sigma^2/2)*t[n]+sigma*W[m,n])
        end
        C[m,n] = exp(-r*t[n])*max(SQ[m,n]-K,0)
        P[m,n] = exp(-r*t[n])*max(K-SQ[m,n],0)
    end
    E[1,n] = mean(S[:,n])
    E[2,n] = S0*exp(mu*t[n])
    Var[1,n] = var(S[:,n],corrected=false)
    Var[2,n] = S0^2*exp(2*mu*t[n])*(exp(sigma^2*t[n])-1)

    Call[1,n] = mean(C[:,n])
    Put[1,n] = mean(P[:,n])

    if n>=2
        d[1,n-1] = (log(S0)-log(K))/(sigma*sqrt(t[n])) + r*sqrt(t[n])/sigma + sigma*sqrt(t[n])/2
        d[2,n-1] = (log(S0)-log(K))/(sigma*sqrt(t[n])) + r*sqrt(t[n])/sigma - sigma*sqrt(t[n])/2
        Call[2,n]= S0*cdf(Normal(),d[1,n-1]) - K*exp(-r*t[n])*cdf(Normal(),d[2,n-1])
        Put[2,n] = K*exp(-r*t[n])*cdf(Normal(),-d[2,n-1]) - S0*cdf(Normal(),-d[1,n-1])
    end
end

plot(t,[E[1,:],E[2,:]], title="True and Estimated Mean of Stock Price",label=["Estimated" "True"],linewidth=2.5)
savefig("True and Estimated Mean of Stock Price.png")
plot(t,[Var[1,:],Var[2,:]], title="True and Estimated Variance of Stock Price",label=["Estimated" "True"],linewidth=2.5)
savefig("True and Estimated Variance of Stock Price.png")
plot(t,[Call[1,:],Call[2,:]], title="True and Estimated Call Price",label=["Estimated" "True"],linewidth=2.5)
savefig("True and Estimated Call Price.png")
plot(t,[Put[1,:],Put[2,:]], title="True and Estimated Put Price",label=["Estimated" "True"],linewidth=2.5)
savefig("True and Estimated Put Price.png")

print("Plots in "*pwd())