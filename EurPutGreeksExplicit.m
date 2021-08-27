%%
% Strike price
K = 100;
% Set the minimal and maximal stock prices
Smin = 0;
Smax = 4*K;
% The number of points in stock direction
N = 1000;
S1 = linspace(Smin,Smax,N+1)';
% The length of the stock price interval
dS = S1(2) - S1(1);
% S stores all the prices except boundary points
S = S1(2:N);

% Maturity of option
T = 1;
% The number of points in time dimension
M = 100000;
tau = linspace(0,T,M+1);
% The length of the time interval
dtau = tau(2) - tau(1);

% Interest rate
r = 0.03;
% Dividend yield
q = 0.05;
% Volatility (Constant)
sigma = 0.2;

% alpha and beta in the discretization scheme
alpha = 0.5*sigma^2*S.^2*dtau/(dS^2);
beta = (r - q)*S*dtau/(2*dS);

% lower, main and upper diagonal of the tridiagonal matrix for the explicit
% finite difference scheme
l = alpha - beta;
d = 1 - r*dtau - 2*alpha;
u = alpha + beta;

%AE = diag(d,0)+diag(l(2:end),-1)+diag(u(1:end-1),1);
% Vector to store the option prices at time k
% Option payoff 
Vold = max(K - S,0);
% Vector to store the option prices at time k+1
Vnew = zeros(N-1,1);

for k=1:M
    % Boundary condition for European put option
    boundary = [l(1)*K*exp(-r*(k-1)*dtau);zeros(N-3,1);u(N-1)*0];
    % Explicit iteration scheme
    for j=1:N-1
        if(j==1)
            Vnew(j) = d(j)*Vold(j) + u(j)*Vold(j+1);
        elseif(j<N-1)
            Vnew(j) = l(j)*Vold(j-1) + d(j)*Vold(j) + u(j)*Vold(j+1);
        else
            Vnew(j) = l(j)*Vold(j-1) + d(j)*Vold(j);
        end
    end
    % Update the vectors from time k to time k+1
    Vold = Vnew + boundary;
end

% Interpolation to find the at-the-money put price when S0 = 100.0
S0 = 100.0;
put_fdm = interp1(S,Vold,S0);

%%
% Compare the prices with the BS formula
[call_bs,put_bs] = blsprice(S,K,r,T,sigma,q);
figure(1)
% Plot the prices from both the explicit method and the BS formula
subplot(2,1,1)
plot(S,Vold,S,put_bs,'LineWidth',2)
title('European Put price, Explicit - BS formula')
xlabel('Stock price')
ylabel('Put price')
legend('Explicit','BS formula','Location','SouthEast')
subplot(2,1,2)
plot(S,Vold - put_bs,'LineWidth',2)
title('Difference of European Put price, Explicit - BS formula')
xlabel('Stock price')
ylabel('Put price')
%%
% Computing the delta of the option using central difference formula:
delta = (Vold(3:1:end) - Vold(1:1:end-2))/(2*dS);
% Computing the delta of the option using Black-Scholes formula:
[call_delta_bs, put_delta_bs] = blsdelta(S(2:end-1),K,r,T,sigma,q);
% Using subplot to plot multiple graphs in one figure
figure(2)
% Plot both delta in one graph on the top pannel
subplot(2,1,1)
plot(S(2:end-1),delta,'LineWidth',2)
hold on;
plot(S(2:end-1),put_delta_bs,'LineWidth',2)
title('European Put delta, Explicit - BS formula')
xlabel('Stock price')
ylabel('Put delta')
legend('Explicit','BS formula','Location','SouthEast')
% Plot the difference between two deltas on the bottom pannel
subplot(2,1,2)
plot(S(2:end-1),delta - put_delta_bs,'LineWidth',2)
title('Differences of European Put delta, Explicit - BS formula')
xlabel('Stock price')
ylabel('Put delta')

%%
% Computing the gamma of the option using central difference formula:
gamma = diff(Vold,2)/(dS*dS);
% Computing the delta of the option using Black-Scholes formula:
gamma_bs = blsgamma(S(2:end-1),K,r,T,sigma,q);
figure(3)
% Plot both delta in one graph on the top pannel
subplot(2,1,1)
plot(S(2:end-1),gamma,'LineWidth',2)
hold on;
plot(S(2:end-1),gamma_bs,'LineWidth',2)
title('European Put gamma, Explicit - BS formula')
xlabel('Stock price')
ylabel('Put gamma')
legend('Explicit','BS formula','Location','SouthEast')
% Plot the difference between two deltas on the bottom pannel
subplot(2,1,2)
plot(S(2:end-1), gamma - gamma_bs,'LineWidth',2)
title('Differences of European Put gamma, Explicit - BS formula')
xlabel('Stock price')
ylabel('Put gamma')