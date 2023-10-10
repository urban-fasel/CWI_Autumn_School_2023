%% Ensemble SINDy 
%% 
% * *Robust model discovery in the high noise and low data limit*
% * *Uncertainty quantification*

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
%% 
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


%%
% library bagging
nEnsemble1P = 0.9; % percentage of full library that is sampled without replacement for library bagging
nEnsemble2 = 100; % number of bootstraps (SINDy models using sampled library terms) in ensemble
ensT = 0.4; % Threshold library term inclusion probabilty: cut library entries that occur less than ensT

% double bagging
nEnsemblesDD = 100; % number of models in ensemble for data bagging after library bagging

% randomly sample library terms without replacement and throw away terms
% with low inclusion probability
nEnsemble1 = round(nEnsemble1P*size(Theta,2));
mOutBS = zeros(nEnsemble1,n,nEnsemble2);
libOutBS = zeros(nEnsemble1,nEnsemble2);
for iii = 1:nEnsemble2
    rs = RandStream('mlfg6331_64','Seed',iii); 
    libOutBS(:,iii) = datasample(rs,1:size(Theta,2),nEnsemble1,'Replace',false)';
    mOutBS(:,:,iii) = sparsifyDynamics(Theta(:,libOutBS(:,iii)),dx,lambda,n);
end

inclProbBS = zeros(size(Theta,2),n);
for iii = 1:nEnsemble2
    for jjj = 1:n
        for kkk = 1:nEnsemble1
            if mOutBS(kkk,jjj,iii) ~= 0
                inclProbBS(libOutBS(kkk,iii),jjj) = inclProbBS(libOutBS(kkk,iii),jjj) + 1;
            end
        end
    end
end
inclProbBS = inclProbBS/nEnsemble2*size(Theta,2)/nEnsemble1;

XiD = zeros(size(Theta,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;
    XiBias = sparsifyDynamics(Theta(:,libEntry),dx(:,iii),lambda,1);
    XiD(libEntry,iii) = XiBias;
end

%%
XiDB = zeros(size(Theta,2),n);
XiDBmed = zeros(size(Theta,2),n);
XiDBs = zeros(size(Theta,2),n);
XiDBeOut = zeros(size(Theta,2),n,nEnsemblesDD);
inclProbDB = zeros(size(Theta,2),n);
for iii = 1:n
    libEntry = inclProbBS(:,iii)>ensT;

    bootstatDD = bootstrp(nEnsemblesDD,@(T,d)sparsifyDynamics(T,d,lambda,1),Theta(:,libEntry),dx(:,iii)); 
    
    XiDBe = [];
    XiDBnz = [];
    for iE = 1:nEnsemblesDD
        XiDBe(:,iE) = reshape(bootstatDD(iE,:),size(Theta(:,libEntry),2),1);
        XiDBnz(:,iE) = XiDBe(:,iE)~=0;
        
        XiDBeOut(libEntry,iii,iE) = XiDBe(:,iE);
    end

    % Thresholded bootstrap aggregating (bagging, from bootstrap aggregating)
    XiDBnzM = mean(XiDBnz,2); % mean of non-zero values in ensemble
    inclProbDB(libEntry,iii) = XiDBnzM;
    XiDBnzM(XiDBnzM<ensembleT) = 0; % threshold: set all parameters that have an inclusion probability below threshold to zero

    XiDBmean = mean(XiDBe,2);
    XiDBmedian = median(XiDBe,2);
    XiDBstd = std(XiDBe')';

    XiDBmean(XiDBnzM==0)=0; 
    XiDBmedian(XiDBnzM==0)=0; 
    XiDBstd(XiDBnzM==0)=0; 
    
    XiDB(libEntry,iii) = XiDBmean;
    XiDBmed(libEntry,iii) = XiDBmedian;
    XiDBs(libEntry,iii) = XiDBstd;
    
end
%%
lib = poolDataLIST({'x','y','z'},XiDBmed,n,polyorder);    
lib(1,:) = [];
skipLastRows = size(XiDBmed,1)-10; % only plot up to quadratic terms
plotUQ(XiDBeOut(1:end-skipLastRows,:,:),true_weights(1:end-skipLastRows,:),XiDBmed,lib)