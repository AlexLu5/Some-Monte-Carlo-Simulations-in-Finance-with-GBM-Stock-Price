S0 = 1;
mu = 0.25;
sigma = 0.5;
r = 0.15;
K = 1;

T = 1;
N = 10000;
M = 10000;
rng("default")

t = linspace(0,T,N*T+1);
Z = normrnd(0,1,M,N*T);
W = zeros(M,N*T+1);
S = zeros(M,N*T+1);
SQ = zeros(M,N*T+1);
C = zeros(M,N*T+1);
P = zeros(M,N*T+1);

E = zeros(2,N*T+1);
Var = zeros(2,N*T+1);
d = zeros(2,N*T);
Call = zeros(2,N*T+1);
Put = zeros(2,N*T+1);

S(:,1) = S0;
SQ(:,1) = S0;
Call(2,1) = max(S0-K,0);
Put(2,1) = max(K-S0,0);

for n = 1:(N*T+1)
    for m = 1:M
        if n>=2
            W(m,n) = W(m,n-1)+1/sqrt(N)*Z(m,n-1);
            S(m,n) = S0*exp((mu-sigma^2/2)*t(n)+sigma*W(m,n));
            SQ(m,n) = S0*exp((r-sigma^2/2)*t(n)+sigma*W(m,n));
        end
        C(m,n) = exp(-r*t(n))*max(SQ(m,n)-K,0);
        P(m,n) = exp(-r*t(n))*max(K-SQ(m,n),0);
    end
    E(1,n) = mean(S(:,n));
    E(2,n) = S0*exp(mu*t(n));
    Var(1,n) = var(S(:,n),1);
    Var(2,n) = S0^2*exp(2*mu*t(n))*(exp(sigma^2*t(n))-1);

    Call(1,n) = mean(C(:,n));
    Put(1,n) = mean(P(:,n));

    if n>=2
        d(1,n-1) = (log(S0)-log(K))/(sigma*sqrt(t(n))) + r*sqrt(t(n))/sigma + sigma*sqrt(t(n))/2;
        d(2,n-1) = (log(S0)-log(K))/(sigma*sqrt(t(n))) + r*sqrt(t(n))/sigma - sigma*sqrt(t(n))/2;
        Call(2,n) = S0*normcdf(d(1,n-1)) - K*exp(-r*t(n))*normcdf(d(2,n-1));
        Put(2,n) = K*exp(-r*t(n))*normcdf(-d(2,n-1)) - S0*normcdf(-d(1,n-1));
    end
end

plot(t,E(1,:),t,E(2,:),"-.",LineWidth=2), xlabel("t"), title("True and Estimated Mean of Stock Price"), legend("Estimated Mean", "True Mean","Location","northwest")
exportgraphics(gcf, "True and Estimated Mean of Stock Price.png")
plot(t,Var(1,:),t,Var(2,:),"-.",LineWidth=2), xlabel("t"), title("True and Estimated Variance of Stock Price"), legend("Estimated Variance", "True Variance","Location","northwest")
exportgraphics(gcf, "True and Estimated Variance of Stock Price.png")
plot(t,Call(1,:),t,Call(2,:),"-.",LineWidth=2), xlabel("t"), title("True and Estimated European Call Price"), legend("Estimated Price", "True Price","Location","northwest")
exportgraphics(gcf, "True and Estimated European Call Price.png")
plot(t,Put(1,:),t,Put(2,:),"-.",LineWidth=2), xlabel("t"), title("True and Estimated European Put Price"), legend("Estimated Price", "True Price","Location","northwest")
exportgraphics(gcf, "True and Estimated European Put Price.png")

disp(strcat("Plots in ",pwd))