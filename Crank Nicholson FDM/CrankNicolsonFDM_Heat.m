% usingCrank Nicolson method to solve the heat equation for u(x,t):
% u_t - u_xx = 0, 
% for 0<=x<=1 and 0<=t<=0.1, with initial condition:
% u(x, 0) = sin(2*pi*x)
% boundary condition
% u(0,t) = 0, u(1,t) = 0.

% Parameters of the above problem:
kappa = 1.0;
a = 0.0;
b = 1.0;
T = 0.1;

% Setup the computational domain:
% The number of time steps
M = 20;
% The time interval
dt = T/M;
% The number of space steps
N = 10; 
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

figure(1)

subplot(2,2,1)
surf(x,t,exact_solution')
xlabel('x')
ylabel('t')
zlabel('u(x,t)')

subplot(2,2,2)
surf(x,t,UCN')
xlabel('x')
ylabel('t')
zlabel('U(x,t)')

%Compute the error between the implicit solution and the exact solution:
error_CN = UCN - exact_solution;
subplot(2,2,3)
surf(x,t,error_CN')
xlabel('x')
ylabel('t')
zlabel('Difference')