function [C_mc, SE_C] = MC_Asian_Call_VR(S0,K,r,D,sigma,T,tm,N_mc)
% % Monte Carlo Simulation to price European Asian call option under GBM
% S0: Stock price at time 0
% K: Strike price of the asian option
% r: risk free interest rate
% D: dividend yield
% sigma: volatility
% T: Maturity of option
% tm, vector of set of dates to compute the average
% N_mc: Number of simulations

% Number of fixed dates from 0=t0<t1<t2<...<tm
m = length(tm)-1;
% Compute each interval of the fixed set of dates
dtm = diff(tm);
% Vector to store sample of the path of stock prices before maturity.
S = zeros(m+1,1);
Sa = S;
% Vector to store sample of discounted call payoff
X = zeros(N_mc,2);
for i=1:N_mc
    S(1) = S0;
    Sa(1) = S0;
    % To simulate and have one sample path S0 = S(1), S(2), S(3), ...
    % S(m+1) and the average is on (S(2), S(3), ... S(m+1))
    Z = randn(m,1);
    for k=1:m
        %Z = randn;
        S(k+1) = S(k)*exp((r-D-0.5*sigma^2)*dtm(k)+sigma*sqrt(dtm(k))*Z(k));
        % Compute the average 
        X(i,1) = ((k-1)*X(i,1)+S(k+1))/(k);
        % Antithetic variate
        Sa(k+1) = Sa(k)*exp((r-D-0.5*sigma^2)*dtm(k)-sigma*sqrt(dtm(k))*(Z(k)));
        X(i,2) = ((k-1)*X(i,2)+Sa(k+1))/(k);
    end
    % Compute the payoff of the asian call option
    X(i,:) = exp(-r*T)*max(X(i,:)-K,0);
end
% The European asian call price is approximated by the sample mean of X
C_mc = mean(mean(X));
% Compute the standard error of the sample
SE_C = sqrt((sum(mean(X,2).^2)-N_mc*C_mc^2)/(N_mc*(N_mc-1)));
% Compute the confidence interval.
%CI_low = C_mc - 1.96*SE_C;
%CI_up = C_mc + 1.96*SE_C;
end

