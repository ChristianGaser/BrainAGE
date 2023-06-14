function [mu, vr, model_gpr] = BA_gpr_core(x_train, y_train, x_test, meanhyp, likhyp)
% [mu, vr, model_gpr] = BA_gpr_core(x_train, y_train, x_test, meanhyp, likhyp)
% Gaussian Process Regression core function
%
% This function implements Gaussian Process Regression (GPR) with a linear
% covariance function using the conjugate gradient method for numerical
% optimization. Given a set of training inputs (x_train) and corresponding
% outputs (y_train), and a set of test inputs (x_test), it predicts the
% mean of the test outputs (mu) using GPR.
%
% INPUTS:
% x_train:   - a matrix of size m x n, where m is the number of training
%              inputs and n is the number of features.
% y_train:   - a vector of length m, corresponding to the training outputs.
% x_test:    - a matrix of size p x n, where p is the number of test inputs
%              and n is the number of features.
% meanhyp:   - mean hyperparameter
% likhyp:    - likelihood hyperparameter for noise variance
%
% OUTPUTS:
% mu:        - a vector of length p, corresponding to the predicted output 
%              mean of the test outputs.
% vr:        - a vector of length p, corresponding to the predicted variances
%              of the test outputs.
% model_gpr  - model structure for predicting new data
%
% The code uses many ideas from the great gpml toolbox by Carl Edward Rasmussen 
% and Hannes Nickisch, but tries to implement only the really necessary steps
% to speed up and save memory.
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if nargin < 3
  error('At least 3 parameters have to be defined. Syntax: mu = BA_gpr_core(x_train, y_train, x_test, meanhyp, likhyp)');
end

% Define covariance function
% Here we use a simple linear covariance function
covfunc = @(x, y) x*y';

% Define prior (constant) mean function
meanfunc = @(x, mn) mn*ones(size(x,1),1);

n_subj = numel(y_train);

% Transpose input if necessary
xs = size(x_train,1);
if xs ~= n_subj
  x_train = x_train';
  x_test = x_test';
end

n_voxel = size(x_train,2);

% Set hyperparameters for mean and likelihood
if nargin < 4
  meanhyp = 100;
end
if nargin < 5
  likhyp  = log10(0.1); % noise variance
end

% Estimate kernel output for training data
K = covfunc(x_train, x_train);  
  
sn2 = exp(2*likhyp);                                % noise variance of likGauss
sW  = sqrt(1/sn2);

% Cholesky factor of B
L = chol(K/sn2 + eye(size(K,1)));
alpha = bsxfun(@times,L\(L'\bsxfun(@times,y_train-meanfunc(x_train, meanhyp),sW)),sW);
clear L

% initialize seed generator
if exist('rng','file') == 2
  rng('default')
  rng(1)
else
  rand('state',1);
end
  
Ks = covfunc(x_train, x_test);  
ms = meanfunc(x_test, meanhyp);                           % evaluate mean vector
mu = ms + Ks'*alpha;                                     % conditional mean fs|f

% prepare predictive output variance
if nargout > 1
  K_post = K + sn2*eye(size(K));
  vr = diag(covfunc(x_test, x_test) - Ks'*((K_post\Ks)));
end

% prepare output of model structure
if nargout > 2
  model_gpr = struct('alpha',alpha,'ms',ms);  
end
