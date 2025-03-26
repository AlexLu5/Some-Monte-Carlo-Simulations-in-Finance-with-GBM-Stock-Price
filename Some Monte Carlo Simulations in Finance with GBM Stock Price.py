import math
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import scipy.stats as ss

S0 = 1
mu = 0.25
sigma = 0.5
r = 0.15
K = 1

T = 1
M = 10000
N = 10000
random.seed(0)

t = np.linspace(0,T,N*T+1)
Z = np.random.normal(0,1,(M,N*T+1))
W = np.empty((M,N*T+1))
S = np.empty((M,N*T+1))
SQ = np.empty((M,N*T+1))
C = np.empty((M,N*T+1))
P = np.empty((M,N*T+1))

E = np.empty((2,N*T+1))
Var = np.empty((2,N*T+1))
d = np.empty((2,N*T+1))
Call = np.empty((2,N*T+1))
Put = np.empty((2,N*T+1))

W[:,0] = 0
S[:,0] = S0
SQ[:,0] = S0
Call[1,0] = max(S0-K,0)
Put[1,0] = max(K-S0,0)

for n in range(0,N*T+1):
    for m in range(0,M):
        if n>=1:
            W[m,n] = W[m,n-1]+1/math.sqrt(N)*Z[m,n-1]
            S[m,n] = S0*math.exp((mu-sigma**2/2)*t[n]+sigma*W[m,n])
            SQ[m, n] = S0*math.exp((r-sigma**2/2)*t[n]+sigma*W[m,n])

        C[m,n] = math.exp(-r*t[n])*max(0,SQ[m,n]-K)
        P[m,n] = math.exp(-r*t[n])*max(0,K-SQ[m,n])

    E[0,n] = np.mean(S[:,n])
    E[1,n] = S0*math.exp(mu*t[n])
    Var[0,n] = np.var(S[:,n])
    Var[1,n] = S0**2*math.exp(2*mu*t[n])*(math.exp(sigma**2*t[n])-1)

    Call[0,n] = np.mean(C[:,n])
    Put[0,n] = np.mean(P[:,n])

    if n>=1:
        d[0,n-1] = (np.log(S0)-np.log(K))/(sigma*math.sqrt(t[n])) + r*math.sqrt(t[n])/sigma + sigma*math.sqrt(t[n])/2
        d[1,n-1] = (np.log(S0)-np.log(K))/(sigma*math.sqrt(t[n])) + r*math.sqrt(t[n])/sigma - sigma*math.sqrt(t[n])/2
        Call[1,n] = S0*ss.norm.cdf(d[0,n-1]) - K*math.exp(-r*t[n])*ss.norm.cdf(d[1,n-1])
        Put[1,n] = K*math.exp(-r*t[n])*ss.norm.cdf(-d[1,n-1]) - S0*ss.norm.cdf(-d[0,n-1])

plt.plot(t,E[0,:],t,E[1,:])
plt.title("True and Estimated Mean of Stock Price")
plt.xlabel("t")
plt.legend(["Estimated","True"])
plt.savefig("True and Estimated Mean of Stock Price")
plt.close()

plt.plot(t,Var[0,:],t,Var[1,:])
plt.title("True and Estimated Variance of Stock Price")
plt.xlabel("t")
plt.legend(["Estimated","True"])
plt.savefig("True and Estimated Variance of Stock Price")
plt.close()

plt.plot(t,Call[0,:],t,Call[1,:])
plt.title("True and Estimated Call Price")
plt.xlabel("t")
plt.legend(["Estimated","True"])
plt.savefig("True and Estimated Call Price")
plt.close()

plt.plot(t,Put[0,:],t,Put[1,:])
plt.title("True and Estimated Put Price")
plt.xlabel("t")
plt.legend(["Estimated","True"])
plt.savefig("True and Estimated Put Price")
plt.close()

print("Plots in "+os.getcwd())