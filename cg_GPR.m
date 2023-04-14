function [mu, vr, model_gpr] = cg_GPR(x_train, y_train, x_test, meanhyp, likhyp, p_dropout, mapping, x_test_org)
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
% p_dropout: - dropout probability to randomly exclude voxels/data points to
%              implement an uncertainty-aware approach using a Monte-Carlo 
%              Dropout during inference. That means that during testing, voxels 
%              are randomly dropped out according to the dropout probabilities.
%              This process is repeated multiple times, and each time, the model 
%              produces a different output. By averaging these outputs, we can 
%              obtain a more robust prediction and estimate the model's 
%              uncertainty in its predictions. A meaningful dropout probability
%              is 0.1, which means that 10% of the data points are excluded.
% mapping    - PCA mapping of dropout is enabled for PCA approach. Here, we have 
%              to use the unmapped data, exclude dropout voxels and apply PCA
%              mapping to the remaining voxels of the test data
% x_test_org - unmapped (original test inputs)
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
  error('At least 3 parameters have to be defined. Syntax: mu = cg_GPR(x_train, y_train, x_test, meanhyp, likhyp)');
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

% Dropout probability
if nargin < 6
  p_dropout  = 0; % no dropout by default
end

% Number of Monte Carlo samples
if p_dropout
  n_mc = 10;
else
  n_mc = 1;
end

% if mapping and original data are defined apply mapping for dropped out data
if nargin == 8
  use_PCA = 1;
else
  use_PCA = 0;
end

% Estimate kernel output for training data
K = covfunc(x_train, x_train);  
  
sn2 = exp(2*likhyp);                                % noise variance of likGauss
sW  = sqrt(1/sn2);

% Cholesky factor of B
L = chol(K/sn2 + eye(size(K,1)));
alpha = bsxfun(@times,L\(L'\bsxfun(@times,y_train-meanfunc(x_train, meanhyp),sW)),sW);
clear L

% Monte Carlo loop
mu_mc = zeros(size(x_test, 1), n_mc);
for i = 1:n_mc

  % initialize seed generator
  if exist('rng','file') == 2
    rng('default')
    rng(i)
  else
    rand('state',i);
  end
    
  % Estimate kernel output for covariance between training and test data
  % using Monte Carlo Dropout approach
  if p_dropout
    
    % for PCA we have to use the unmapped data, exclude dropout voxels and
    % apply PCA mapping to the remaining voxels of the test data
    if use_PCA

      % exclude random voxels according to dropout probability
      ind = randperm(size(x_test_org,2)) > p_dropout*size(x_test_org,2);

      % apply PCA mapping to remaining voxels of test data
      mapped_test = (x_test_org(:,ind) - repmat(mapping.mean(ind), [size(x_test_org, 1) 1])) * mapping.M(ind,:);
      mapped_test  = (mapped_test - mapping.mn)/(mapping.mx - mapping.mn);
      Ks = covfunc(x_train, mapped_test);
    else
      % exclude random voxels according to dropout probability
      ind = randperm(n_voxel) > p_dropout*n_voxel;
      Ks = covfunc(x_train(:,ind), x_test(:,ind));
    end
  else
    Ks = covfunc(x_train, x_test);  % no indexing is much faster without dropout
  end
  ms = meanfunc(x_test, meanhyp);                         % evaluate mean vector
  mu_mc(:, i) = ms + Ks'*alpha;                          % conditional mean fs|f
end

% Compute final prediction by taking the mean of the Monte Carlo samples
mu = mean(mu_mc, 2);

% prepare predictive output variance
if nargout > 1
  K_post = K + sn2*eye(size(K));
  vr = diag(covfunc(x_test, x_test) - Ks'*((K_post\Ks)));
end

% prepare output of model structure
if nargout > 2
  model_gpr = struct('alpha',alpha,'ms',ms);  
end
