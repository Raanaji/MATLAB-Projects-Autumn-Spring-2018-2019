function [C_mc, SE_C] = MC_Binary_Put_VR(S0,K,r,D,sigma,T,N_mc)
% % Monte Carlo Simulation to price European binary put option under GBM
% S0: Stock price at time 0
% K: Strike price of the asian option
% r: risk free interest rate
% D: dividend yield
% sigma: volatility
% T: Maturity of option
% N_mc: Number of simulations

% Vector to store sample of stock prices at maturity
S = zeros(N_mc,1);
% Vector to store sample of discounted call payoff
X = zeros(N_mc,2);
%Y = X;
for i=1:N_mc
    Z = randn;
    S(i) = S0*exp((r-D-0.5*sigma^2)*T+sigma*sqrt(T)*Z);
    % Payoff of the binary put option
    X(i,1) = exp(-r*T)*(S(i)<=K);
    % Antithetic variates
    S(i) = S0*exp((r-D-0.5*sigma^2)*T+sigma*sqrt(T)*(-Z));
    % Payoff of the binary put option
    %Y(i) = exp(-r*T)*(S(i)<=K);
    X(i,2) = exp(-r*T)*(S(i)<=K);
    %X(i) = 0.5*(X(i) + Y(i));
end
% The European binary put price is approximated by the sample mean of X
C_mc = mean(mean(X));
% Compute the standard error of the sample
SE_C = sqrt((sum(mean(X,2).^2)-N_mc*C_mc^2)/(N_mc*(N_mc-1)));
% Compute the confidence interval.
%CI_low = C_mc - 1.96*SE_C;
%CI_up = C_mc + 1.96*SE_C;
end