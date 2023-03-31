function mu = cg_GPR(x_train, y_train, x_test, meanhyp, likhyp, p_dropout)
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
%
% OUTPUTS:
% mu:        - a vector of length p, corresponding to the predicted mean of the
%              test outputs.
%
% The code uses many ideas from the great gpml toolbox by Carl Edward Rasmussen 
% and Hannes Nickisch, but tries to implement only the really necessary steps
% to speed up and save memory.

if nargin < 3
  error('At least 3 parameters have to be defined. Syntax: mu = cg_GPR(x_train, y_train, x_test, meanhyp, likhyp)');
end

n_subj = numel(y_train);

% Transpose input if necessary
xs = size(x_train,1);
if xs ~= n_subj
  x_train = x_train';
  x_test = x_test';
end

n_voxel = size(x_train,2);

% Define covariance function
% Here we use a simple linear covariance function
covfunc = @(x, y) x*y';

% Define prior (constant) mean function
meanfunc = @(x, mn) mn*ones(size(x,1),1);

% Set hyperparameters for mean and likelihood
if nargin < 4
  meanhyp = 100;
end
if nargin < 5
  likhyp  = log10(0.1); % noise variance
end

% Dropout probability
if nargin < 6
  p_dropout  = 0; % no deropout by default
end

% Number of Monte Carlo samples
if p_dropout
  n_mc = 10;
else
  n_mc = 1;
end

% Estimate kernel output for training data
K = covfunc(x_train, x_train);  
  
sn2 = exp(2*likhyp);                                % noise variance of likGauss
W   = ones(size(K,1),1)/sn2;

% Go through Monte Carlo
mu_mc = zeros(size(x_test, 1), n_mc);
for i = 1:n_mc

  % initialize seed generator
  if exist('rng','file') == 2
    rng('default')
    rng(i)
  else
    rand('state',i);
  end
    
  % switch between Cholesky and LU decomposition mode
  isWneg = any(W<0); n = numel(W);
  if isWneg 
    A = eye(n) + bsxfun(@times,K,W');                     % asymmetric A = I+K*W
    [L,U,P] = lu(A); u = diag(U);         % compute LU decomposition, A = P'*L*U
    detP = 1;               % compute sign (and det) of the permutation matrix P
    p = P*(1:n)';
    for j=1:n                                                     % swap entries
        if j~=p(j), detP = -detP; k = find(p==j); p([j,k]) = p([k,j]); end
    end
  else                                                 % symmetric B = I+sW*K*sW
    sW = sqrt(W); L = chol(eye(n)+sW*sW'.*K);             % Cholesky factor of B
    solveKiW = @(r) bsxfun(@times,L\(L'\bsxfun(@times,r,sW)),sW);
  end

  alpha = solveKiW(y_train-meanfunc(x_train, meanhyp));

  % Estimate kernel output for covariance between training and test data
  if p_dropout
    % exclude random voxels according to dropout probability
    ind = randperm(n_voxel) > p_dropout*n_voxel;
    Ks = covfunc(x_train(:,ind), x_test(:,ind));
  else
    Ks = covfunc(x_train, x_test);  % no indexing is much faster without dropout
  end
  ms = meanfunc(x_test, meanhyp);                         % evaluate mean vector
  mu_mc(:, i) = ms + Ks'*alpha;                          % conditional mean fs|f
end

% Compute final prediction by taking the mean of the Monte Carlo samples
mu = mean(mu_mc, 2);

