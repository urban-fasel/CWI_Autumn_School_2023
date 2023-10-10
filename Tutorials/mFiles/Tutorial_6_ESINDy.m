%% Ensemble SINDy 
% *Robust model discovery* in the high noise and low data limit
% 
% *Uncertainty quantification*: 
%% 
% * Model structure: coefficient inclusion probability
% * Model parameter uncertainty
%% 
% Paper: Fasel et al: <https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2021.0904 
% https://royalsocietypublishing.org/doi/abs/10.1098/rspa.2021.0904> 
% 
% 

addpath(genpath(pwd))
%% Lorenz system

% Generate data
param = [10; 28; 8/3]; % Lorenz system parameters (chaotic)
dt = 0.01; % time step
tFinal = 10; % final time
tspan = dt:dt:tFinal;
n = 3; % number of states
x0 = [-8; 8; 27];  % Initial condition
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,3));
[tClean,xClean]=ode45(@(t,x) lorenz(t,x,param),tspan,x0,options);

% true weights to compare results
true_weights = zeros(20,n);
true_weights(2,1) = -10;
true_weights(3,1) = 10;
true_weights(2,2) = 28;
true_weights(3,2) = -1;
true_weights(7,2) = -1;
true_weights(4,3) = -8/3;
true_weights(6,3) = 1;

%% E-SINDy vs SINDy

% ensemble hyperparameters: data bootstraps (vs library bagging)
nEnsembles = 100; % number of bootstraps (SINDy models using sampled data) in ensemble
ensembleT = 0.6; % Threshold model coefficient inclusion probability: set ensemble SINDy model coefficient to zero if inclusion probability is below ensembleT

% Add Gaussian white noise
sig_nn = 0:0.01:0.05; % noise level
n_noise = length(sig_nn);
m_noise = 5; % number of noise realizations

% output
nWrongTerms = zeros(n_noise,m_noise);
modelError = zeros(n_noise,m_noise);
success = zeros(n_noise,m_noise);
nWrongTermsE = zeros(n_noise,m_noise);
modelErrorE = zeros(n_noise,m_noise);
successE = zeros(n_noise,m_noise);

for mm = 1:m_noise
    for nn = 1:n_noise
        
        % Add noise
        rng(mm) % specify seed for the random number generator for reproducibility
        x = xClean + sig_nn(nn)*std(xClean(:))*randn(size(xClean));
    
        % Compute Derivative: finite difference
        dx = (1/(12*dt))*(-x(5:end,:)+8*x(4:end-1,:)-8*x(2:end-3,:)+x(1:end-4,:)); % fourth order central difference
        x = x(3:end-2,:); % cut tails
        
        % Pool Data (i.e., build library of nonlinear time series)\
        polyorder = 3; % polynomials up to order 3
        Theta = poolData(x,n,polyorder);
        
        % SINDy
        lambda = 0.2; % lambda is our sparsification knob.
        Xi = sparsifyDynamics(Theta,dx,lambda,n); % identify model coefficients
        
        % E-SINDy: identify SINDy models on bootstraps of data
        [bootstat,bootstatn] = bootstrp(nEnsembles,@(T,d)sparsifyDynamics(T,d,lambda,n),Theta,dx); 
        XiE = [];
        XiEnz = [];
        for iE = 1:nEnsembles
            XiE(:,:,iE) = reshape(bootstat(iE,:),size(Theta,2),n);
            XiEnz(:,:,iE) = XiE(:,:,iE)~=0;
        end
    
        % Aggregate SINDy models: inclusion probability
        XiEnzM = mean(XiEnz,3); % mean of non-zero values in ensemble
        XiEnzM(XiEnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero
        
        % bragging SINDy: median of model coefficients
        XiMedian = median(XiE,3);
        XiMedian(XiEnzM==0)=0; 
        
        % number of wrong terms and model coefficient error
        nWrongTerms(nn,mm) = sum(sum(abs((true_weights~=0) - (Xi~=0))));
        modelError(nn,mm) = norm(Xi-true_weights)/norm(true_weights);    
        nWrongTermsE(nn,mm) = sum(sum(abs((true_weights~=0) - (XiMedian~=0))));
        modelErrorE(nn,mm) = norm(XiMedian-true_weights)/norm(true_weights);
    end
end

% average over noise realizations
nWrongTerms = mean(nWrongTerms,2);
modelError = mean(modelError,2);
nWrongTermsE = mean(nWrongTermsE,2);
modelErrorE = mean(modelErrorE,2);

plotESINDy(sig_nn,nWrongTerms,nWrongTermsE,modelError,modelErrorE)

% Plot E-SINDy model coefficient distributions
% Legend plot
%% 
% * coloured subplot boxes:    true model (true active terms)
% * magenta cross:                                true model coefficients 
% * coloured horizontal line:    identified model (identified active terms)

lib = poolDataLIST({'x','y','z'},XiMedian,n,polyorder);    
lib(1,:) = [];
skipLastRows = size(XiMedian,1)-10; % only plot up to quadratic terms
plotUQ(XiE(1:end-skipLastRows,:,:),true_weights(1:end-skipLastRows,:),XiMedian,lib)

%% Evaluate the median models from the two clusters in the bimodal distribution

% split the two models from two clusters in the bimodal distribution
mm = [squeeze(XiE(2,2,:)),squeeze(XiE(3,2,:))]; % two coefficients that show bimodal distribution
GMModel = fitgmdist(mm,2); % fit a Gaussian mixture distribution to data
clusterXi = cluster(GMModel,mm); % Cluster index 

% bragging: take the median of the two clusters
Xi1 = median(XiE(:,:,clusterXi==1),3);
Xi2 = median(XiE(:,:,clusterXi==2),3);

% test and plot model 1
paramSINDy.Xi = Xi1;
paramSINDy.polyorder = polyorder;
[tSINDy,xSINDy]=ode45(@(t,x) SINDyODE(t,x,paramSINDy),tspan,x0,options);
plotSINDy(tClean,xClean,tSINDy(3:end-2,:),xSINDy(3:end-2,:))

% test and plot model 2
paramSINDy.Xi = Xi2;
paramSINDy.polyorder = polyorder;
[tSINDy,xSINDy]=ode45(@(t,x) SINDyODE(t,x,paramSINDy),tspan,x0,options);
plotSINDy(tClean,xClean,tSINDy(3:end-2,:),xSINDy(3:end-2,:))