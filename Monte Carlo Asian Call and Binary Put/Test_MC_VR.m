% Paramters of GBM
r = 0.05;
D = 0.0;
sigma = 0.2;

% Specifications of the option
S0 = 50;
K = 40;
T = 1;
M = 12;
dt = T/M;
tm = 0:dt:T;

% Number of Simulations
N_mc = 10000;

% Price of a Binary put option
[Price_BinaryPut, SE_BinaryPut] = MC_Binary_Put(S0,K,r,D,sigma,T,N_mc);

% Antithetic Variates
[Price_BinaryPut_VR, SE_BinaryPut_VR] = MC_Binary_Put_VR(S0,K,r,D,sigma,T,N_mc);

% Control Variates
[Price_BinaryPut_CV, SE_BinaryPut_CV] = MC_Binary_Put_CV(S0,K,r,D,sigma,T,N_mc);

% Price of an Asian call option
[Price_AsianCall, SE_AsianCall] = MC_Asian_Call(S0,K,r,D,sigma,T,tm,N_mc);

% Antithetic Variates
[Price_AsianCall_VR, SE_AsianCall_VR] = MC_Asian_Call_VR(S0,K,r,D,sigma,T,tm,N_mc);

% Analytic price of the Binary put option:
dm = (log(S0/K)+(r-D-0.5*sigma^2)*T)/(sigma*sqrt(T));
Price_BinaryPut_Formula = exp(-r*T)*(1 - cdf('normal',dm,0,1));


%% 
% Test the basket exchange option
d = 1;
r = 0.05;
T = 1;
c = ones(d,1)/d;
rho = -0.25;
sigma = [0.2;0.3];
Sig = [1 rho;rho 1];
S0 = [50;50];
%sigma = [0.2;0.2;0.3;0.3];
%Sig = [1 0.1 0.2 0.2; 0.1 1 0.3 -0.3; 0.2 0.3 1 0.5;0.2 -0.3 0.5 1];
%S0 = [50*ones(d,1);45];
n = 1000000;

[price, se] = MC_Exchange_Option(r,T,sigma,c,Sig,S0,n)

%%
% Analytic price of exchange option under GBM
sig = sqrt(sum(sigma.^2)-2*rho*sigma(1)*sigma(2));
dp = (log(S0(1)/S0(2)) + 0.5*sig^2*T)/(sig*sqrt(T));
dm = (log(S0(1)/S0(2)) - 0.5*sig^2*T)/(sig*sqrt(T));
M = S0(1)*cdf('normal',dp,0,1) - S0(2)*cdf('normal',dm,0,1)