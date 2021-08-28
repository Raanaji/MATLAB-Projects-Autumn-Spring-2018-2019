function [price, se] = MC_Exchange_Option(r,T,sigma,c,Sig,S0,n)
% A function to compute the prices and standard error of a basket exchange
% option
% r - risk free rate
% T - Maturity
% sigma - vector of volatilities of each asset
% c - weights on the baskets
% Sig - covariance matrix
% S0 - vector of initial stock prices
% n = number of simulations

% d - dimension of the problem
d = length(c);

% Cholesky decomposition of the covariance matrix Sig:
A = chol(Sig,'lower'); 

% ST - a vector of (d+1)*1 the sample prices at maturity
ST = zeros(d+1,1);

% X is a vector of size n to store the payoff for each sample.
X = zeros(n,2);

% Starting the simulation of the stock prices
for i=1:n
    Z = A*randn(d+1,1);
    %for j=1:d+1
    %    ST(j) = S0(j)*exp((r-0.5*sigma(j)^2)*T+sigma(j)*sqrt(T)*Z(j));
    %end
    % Stock prices at maturity
    ST = S0.*exp((r-0.5*sigma.^2)*T+sigma.*sqrt(T).*Z);
    % Payoff at maturity
    X(i,1) = exp(-r*T)*max(sum(c.*ST(1:d)) - ST(d+1),0);
    % Antithetic variates
    ST = S0.*exp((r-0.5*sigma.^2)*T+sigma.*sqrt(T).*(-Z));
    X(i,2) = exp(-r*T)*max(sum(c.*ST(1:d)) - ST(d+1),0);
end
% The exchange option price is approximated by the sample mean of X
price = mean(mean(X));
% Compute the standard error of the sample
se = sqrt((sum(mean(X,2).^2)-n*price^2)/(n*(n-1)));
end