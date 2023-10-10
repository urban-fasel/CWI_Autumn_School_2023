%% SINDy PDE find: Burgers. equation
% Burgers. equation: $\partial_t \;u+u\partial_x \;u=\nu \;/\;\pi \;\;\partial_x^2 
% \;u\;$
% 
% Burgers’ equation is derived from the Navier Stokes equations for the velocity 
% field by dropping the pressure gradient term. 
%% 
% * Despite it’s relation to the much more complicated Navier Stokes equations, 
% Burgers’ equation does not exhibit turbulent behavior. 
% * 'Burgers_data': solution to Burgers’ equation with a Gaussian initial condition, 
% propagating into a traveling wave. 
% * Unlike many solutions to the inviscid Burgers’ equation (ut + uux = 0), 
% the dissipative term uxx prevents a shock from forming. 
%% 
% Code and data modified from Rudy et al (2017) and Reinbold et al (2020).
%% 
% * Paper Rudy: <https://www.science.org/doi/10.1126/sciadv.1602614 https://www.science.org/doi/10.1126/sciadv.1602614> 
% * Github Rudy: <https://github.com/snagcliffs/PDE-FIND https://github.com/snagcliffs/PDE-FIND> 
% * Paper Reinbold: <https://journals.aps.org/pre/pdf/10.1103/PhysRevE.101.010203 
% https://journals.aps.org/pre/pdf/10.1103/PhysRevE.101.010203> 
% * Github Reinbold: <https://github.com/pakreinbold/PDE_Discovery_Weak_Formulation 
% https://github.com/pakreinbold/PDE_Discovery_Weak_Formulation>   
%% Burgers' data

addpath(genpath(pwd))

% load data: Rudy 2017, coefficient nu/pi = 0.1
traj = matfile('Burgers_data.mat'); 

% grid densities
dt = traj.dt;
dx = traj.dx; 

% add noise
sig = 0.0;
% sig = 0.01;

Ua = traj.uu; 
Ua = Ua + sig*std(Ua(:))*randn(size(Ua)); % Gaussian

% plot data
plotPDEdata(traj,Ua,'Burgers equation')

%% 
%% Identify PDE with standard SINDy

% compute derivatives using finite difference
ut=[];
for xpos=1:size(Ua,1)
    ut(xpos,:)=FiniteDiff(real(Ua(xpos,:)),dt,1);
end

ux=[];uxx=[];uxxx=[];uxxxx=[];
for kk=1:size(Ua,2)
    ux(:,kk)=FiniteDiff(real(Ua(:,kk)),dx,1);
    uxx(:,kk)=FiniteDiff(real(Ua(:,kk)),dx,2);
    uxxx(:,kk)=FiniteDiff(real(Ua(:,kk)),dx,3);
    uxxxx(:,kk)=FiniteDiff(real(Ua(:,kk)),dx,4);
end

%%
%% time derivatives: left hand side
LHS = ut(:);

%% Library Theta: right hand side
const = ones(traj.Nt,traj.Nx);
u = Ua;
u2 = Ua.^2; 
u3 = Ua.^3;
% first order derivative ux 
uux = Ua.*ux; % advection term
u2ux = Ua.^2.*ux;
u3ux = Ua.^3.*ux;
% second order derivative uxx: Laplacian 
uuxx = Ua.*uxx;
u2uxx = Ua.^2.*uxx;
u3uxx = Ua.^3.*uxx;
% third order derivative uxxx
uuxxx = Ua.*uxxx;
u2uxxx = Ua.^2.*uxxx;
u3uxxx = Ua.^3.*uxxx;
% fourth order derivative: Biharmonic 
uuxxxx = Ua.*uxxxx;
u2uxxxx = Ua.^2.*uxxxx;
u3uxxxx = Ua.^3.*uxxxx;

Theta = [const(:), u(:), u2(:), u3(:), ux(:), uux(:), u2ux(:), u3ux(:), uxx(:), uuxx(:), u2uxx(:), u3uxx(:), uxxx(:), uuxxx(:), u2uxxx(:), u3uxxx(:), uxxxx(:), uuxxxx(:), u2uxxxx(:), u3uxxxx(:),];

%%
% SINDy using sequentially trhesholding least squares
n = 1; % number of states
lambda = 0.05; % sparsification knob
Xi = sparsifyDynamics(Theta,LHS,lambda,n); % SINDy

% print coefficients
ThetaLIST = {'1','u','u^2','u^3','u_{x}','uu_{x}','u^2u_{x}','u^3u_{x}','u_{xx}','uu_{xx}','u^2u_{xx}','u^3u_{xx}','u_{xxx}','uu_{xxx}','u^2u_{xxx}','u^3u_{xxx}','u_{xxxx}','uu_{xxxx}','u^2u_{xxxx}','u^3u_{xxxx}'}';
Xi_SINDy = [ThetaLIST, num2cell(Xi)]

% error
c = [-1; 0.1]; % KS true model coefficients
XiTrue = Xi([6 9]); % KS true model structure
error_SINDy = abs(XiTrue-c)./c % model coefficient error
nonZeroTermsSINDy = nnz(Xi) % number of nonzeros
  
%%