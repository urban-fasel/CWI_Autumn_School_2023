%% "Vanilla" SINDy limitations / challenges 
%% 
% * *Data*: noise level, amount of data, sampling rate
% * *Model selection*: hyperparameter lambda
% * *Library*: polynomials, trigonometric functions
%% 
%% 1. Data: noise level, amount of data, sampling rate

param = [10; 28; 8/3]; % Lorenz system parameters (chaotic)
n = 3; % number of states
x0 = [-8; 8; 27];  % Initial condition
dt = 0.001; % time step
tFinal = 20; % final time
tspan = dt:dt:tFinal;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[t,x]=ode45(@(t,x) lorenz(t,x,param),tspan,x0,options);


% Add Gaussian white noise
rng(1)
sig =0.01; % noise level
x = x + sig*std(x(:))*randn(size(x));


% Compute Derivative: finite difference
dx = (1/(12*dt))*(-x(5:end,:)+8*x(4:end-1,:)-8*x(2:end-3,:)+x(1:end-4,:)); % fourth order central difference
x = x(3:end-2,:); % cut tails


% Pool Data (i.e., build library of nonlinear time series)
polyorder = 3; % polynomials up to order 3
Theta = poolData(x,n,polyorder);
m = size(Theta,2); % size of library


% Compute Sparse regression: sequential least squares
lambda = 0.1;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n); % identify model coefficients
disp(poolDataLIST({'x','y','z'},Xi,n,polyorder)) % display SINDy model

%% 
%% 2. Model selection: hyperparameter$\lambda$ 
%% Sweep over lambda, starting from LS solution ($\lambda =0$)

lambdaP = -4; 
lambda = 10^lambdaP; % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n); % identify model coefficients
disp(poolDataLIST({'x','y','z'},Xi,n,polyorder)) % display SINDy model

%% 
%% 3. Library: polynomials, trigonometric functions, ... 
%% Generate Data - inverted pendulum

paramP.m = 0.1; % mass     
paramP.l = 0.5; % arm length
paramP.k = 0.03; % friction coefficient
paramP.g = 9.8; % gravity
n = 2; % number of states
dt = 0.001; % time step
tspan = dt:dt:20;
forcing = @(x,t) 0; % pendulum forcing
x0=[0.1;0];  % Initial condition: close to upright position
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
[t,x]=ode45(@(t,x) pendulumODE(t,x,forcing(x,t),paramP),tspan,x0,options);

% Add Gaussian white noise
rng(1)
sig = 0.0;  
x = x + sig*std(x(:))*randn(size(x));

plotPendulum(t,x) % plot data

%% Compute Derivative: finite difference

dx = (1/(12*dt))*(-x(5:end,:)+8*x(4:end-1,:)-8*x(2:end-3,:)+x(1:end-4,:)); % fourth order central difference
x = x(3:end-2,:); % cut tails
t = t(3:end-2,:);

%% Pool Data (i.e., build library of nonlinear time series)

polyorder = 3; % polynomials up to order 3
Theta = poolData(x,n,polyorder);

%% Compute Sparse regression: sequential least squares

lambda = 0.025;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n); % identify model coefficients
disp(poolDataLIST({'x','xD'},Xi,n,polyorder)) % display SINDy model

%%  Use SINDy model for prediction

paramSINDy.Xi = Xi;
paramSINDy.polyorder = polyorder;
[tSINDy,xSINDy]=ode45(@(t,x) SINDyODE(t,x,paramSINDy),tspan,x0,options);

plotSINDyPendulum(t,x,tSINDy(3:end-2,:),xSINDy(3:end-2,:)) 

%% 
%% Use trigonometric functions in library

% build library
polyorder = 3; % polynomials up to order 3
usesine = 1;
sineorder = 2;
Theta = poolDataSine(x,n,polyorder,usesine,sineorder);
m = size(Theta,2); % size of library

%%
% identify SINDy model
lambda = 0.25;      % lambda is our sparsification knob.
Xi = sparsifyDynamics(Theta,dx,lambda,n); % identify model coefficients
disp(poolDataLISTsine({'x','xD'},Xi,n,polyorder,usesine,sineorder)) % display SINDy model

%%
% SINDy model prediction
paramSINDy.Xi = Xi;
paramSINDy.polyorder = polyorder;
paramSINDy.usesine = usesine;
paramSINDy.sineorder = sineorder;
[tSINDy,xSINDy]=ode45(@(t,x) SINDyODEsine(t,x,paramSINDy),tspan,x0,options);

plotSINDyPendulum(t,x,tSINDy(3:end-2,:),xSINDy(3:end-2,:))