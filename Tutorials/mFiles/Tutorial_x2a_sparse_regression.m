%% OLS, ridge regression, LASSO
%% 
% * Ordinary Least Squares:    $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} 
% \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 \;\right.$
% * Ridge Regression:             $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} 
% \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 +\textrm{𝜆}\right\|\textrm{𝐱}{\left\|\right.}_2^{2\;}$
% * LASSO Regression:           $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} 
% \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 +\textrm{𝜆}\right\|\textrm{𝐱}{\left\|\right.}_1$
% * Elastic Net Regression:      $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} 
% \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 +\lambda_1 \right\|\textrm{𝐱}{\left\|\right.}_1 
% +\lambda_2 \left\|\textrm{𝐱}{\left\|\right.}_2^{2\;} \right.$
%% Generate data

addpath(genpath(pwd))

m = 100;    % number of measurements (time steps)
D = 10;     % number of features (library terms)

% Matrix of possible predictors (library in SINDy)
rng(1) % specify seed random number generator
A = randn(m,D); % normal distribution with mean 0 and standard deviation 1                 

% Three nonzero predictors
x = [0; 1; 1; 0; 0; 0; -1; 0; 0; 0];    

% level of collinearity: near perfect multicollinearity for e.g. 0.001
A(:,2) = A(:,3) + 0.001*randn(m,1);    

% Observations (with noise)
b = A*x + 1.0*randn(m,1);                 

% 
% Ordinary Least Squares: $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 \;\right.$ 

xHat_OLS = (A'*A)^-1*A'*b
% 
% Ridge Regression: $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 +\textrm{𝜆}\right\|\textrm{𝐱}{\left\|\right.}_2^{2\;}$

lambda = 0.036; % regularization param: 0 -> OLS
xHat_ridge = (A'*A+lambda*eye(D))^-1*A'*b % usually, first normalize columns of A to mean 0 and std 1
% 
% LASSO Regression: $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 +\textrm{𝜆}\right\|\textrm{𝐱}{\left\|\right.}_1$

% 10-fold cross-validation, coordinate descent algorithm 
[XL1,FitInfo] = lasso(A,b,'CV',10);
xHat_LASSO = XL1(:,FitInfo.Index1SE)
%%
% de-biasing: OLS on selected features 
xHat_LASSO_DeBiased = zeros(D,1);
xHat_LASSO_DeBiased(abs(xHat_LASSO)>0) = A(:,abs(xHat_LASSO)>0)\b
%%
% de-biasing 2: ridge regression on selected features 
xHat_LASSO_Ridge = zeros(D,1);
lambda = 0.025;
xHat_LASSO_Ridge(abs(xHat_LASSO)>0) = (A(:,abs(xHat_LASSO)>0)'*A(:,abs(xHat_LASSO)>0)+lambda*eye(sum(abs(xHat_LASSO)>0)))^-1*A(:,abs(xHat_LASSO)>0)'*b
%% 
% Elastic Net Regression: $\textrm{𝐱}̂={\textrm{argmin}}_{\mathbf{x}} \left\|\textrm{𝐛}-\textrm{𝐀𝐱}{\left\|\right.}_2^2 +\lambda_1 \right\|\textrm{𝐱}{\left\|\right.}_1 +\lambda_2 \left\|\textrm{𝐱}{\left\|\right.}_2^{2\;} \right.$

alpha = 0.5; % alpha = 1 -> LASSO, alpha close to 0 -> ridge regression
[XL2,FitInfoEN] = lasso(A,b,'CV',10,'Alpha',0.5);

xHat_EN = XL2(:,FitInfoEN.Index1SE)