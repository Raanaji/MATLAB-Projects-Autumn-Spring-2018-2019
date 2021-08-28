%%
% Code to test Crank Nicolson method
kappa = 1.0;
a = 0.0;
b = 1.0;
T = 0.1;
M = 20;
N = 10;
UI = CrankNicolsonFDMHeat(kappa,a,b,T,N,M);
surf(UI')
%%
% Test the order of convergence in space
% Let the number of points big enough in time direction
M = 5000;
N1 = 10;
[UI, error_space1] = CrankNicolsonFDMHeat(kappa,a,b,T,N1,M);
N2 = 2*N1;
[UI, error_space2] = CrankNicolsonFDMHeat(kappa,a,b,T,N2,M);
% order p
p_space = log2(error_space1/error_space2)
%%
% Test the order of convergence in time
% Let the number of points big enough in space direction
N = 5000;
M1 = 10;
[UI, error_time1] = CrankNicolsonFDMHeat(kappa,a,b,T,N,M1);
M2 = 2*M1;
[UI, error_time2] = CrankNicolsonFDMHeat(kappa,a,b,T,N,M2);
p_time = log2(error_time1/error_time2)
