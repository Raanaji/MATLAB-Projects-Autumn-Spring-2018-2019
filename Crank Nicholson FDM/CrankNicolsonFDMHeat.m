% using Crank Nicolson method to solve the heat equation for u(x,t):
% u_t - u_xx = 0, 
% for 0<=x<=1 and 0<=t<=0.1, with initial condition:
% u(x, 0) = sin(2*pi*x)
% boundary condition
% u(0,t) = 0, u(1,t) = 0.

function [UCN, error] = CrankNicolsonFDMHeat(kappa,a,b,T,N,M)
% Setup the computational domain:
% The number of time steps : M
% The time interval
dt = T/M;
% The number of space steps : N 
% The space interval
dx = (b-a)/N;
% parameter rho
rho = kappa*dt/dx^2;

% The grid in space
x = a+dx:dx:b-dx;
% The grid in time
t = 0:dt:T;

% Tri-diagonal Matrix for the implicit scheme AI
%AI = diag((1+2*rho)*ones(N-1,1)) + diag(-rho*ones(N-2,1),1)...
%    + diag(-rho*ones(N-2,1),-1);

% Instead of Creating a tri-diagonal matrix, we only need to save three
% vectors for main, upper and lower diagnoal of the matrix:
% diagnoal 
dI = (1+2*rho)*ones(N-1,1);
% upper diagnoal
uI = -rho*ones(N-1,1);
% lower diagnoal
lI = -rho*ones(N-1,1);

dE = (1-2*rho)*ones(N-1,1);
uE = rho*ones(N-1,1);
lE = rho*ones(N-1,1);

% Matrix to store the vector of solution on grid points at each time step
UCN = zeros(N-1,M+1);
% Initial condition
UCN(:,1) = sin(2*pi.*x');

% Iteration in time direction to compute all solutions:
for k=1:M
    % Compute the explicit side of the equation
    rhs = tridiag_prod(dE+1,uE,lE,UCN(:,k));
    % Using tri-diaganol solver to find the value function at k+1
    UCN(:,k+1) = tridiag(dI+1,uI,lI,rhs);
end

exact_solution = zeros(length(x),length(t));
for j=1:length(x)
    exact_solution(j,:) = exp(-4*pi^2*t)*sin(2*pi*(x(j)));
end

%Compute the error between the implicit solution and the exact solution:
error_CN = exact_solution - UCN;
error = sqrt(mean(error_CN(:,M+1).^2));
end