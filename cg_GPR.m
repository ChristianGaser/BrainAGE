function mu = cg_GPR(x_train, y_train, x_test, meanhyp, likhyp)
% This function implements Gaussian Process Regression (GPR) with a linear
% covariance function using the conjugate gradient method for numerical
% optimization. Given a set of training inputs (x_train) and corresponding
% outputs (y_train), and a set of test inputs (x_test), it predicts the
% mean of the test outputs (mu) using GPR.
%
% INPUTS:
% - x_train: a matrix of size m x n, where m is the number of training
%            inputs and n is the number of features.
% - y_train: a vector of length m, corresponding to the training outputs.
% - x_test:  a matrix of size p x n, where p is the number of test inputs
%            and n is the number of features.
% - meanhyp: mean hyperparameter
% - likhyp:  likelihood hyperparameter for noise variance
%
% OUTPUTS:
% - mu:     a vector of length p, corresponding to the predicted mean of the
%           test outputs.
%
% The code uses many ideas from the great gpml toolbox by Carl Edward Rasmussen 
% and Hannes Nickisch, but tries to implement only the really necessary steps
% to speed up and save memory.

if nargin < 3
  error('At least 3 parameters have to be defined. Syntax: mu = cg_GPR(x_train, y_train, x_test, meanhyp, likhyp)');
end

n_subj = numel(y_train);

% if neede transpose input
xs = size(x_train,1);
if xs ~= n_subj
  x_train = x_train';
  x_test = x_test';
end

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

% Estimate kernel output for training data
K = covfunc(x_train, x_train);

sn2 = exp(2*likhyp);                               % noise variance of likGauss
W   = ones(size(K,1),1)/sn2;

% switch between Cholesky and LU decomposition mode
isWneg = any(W<0); n = numel(W);
if isWneg
  A = eye(n) + bsxfun(@times,K,W');                     % asymmetric A = I+K*W
  [L,U,P] = lu(A); u = diag(U);         % compute LU decomposition, A = P'*L*U
  detP = 1;               % compute sign (and det) of the permutation matrix P
  p = P*(1:n)';
  for i=1:n                                                     % swap entries
    if i~=p(i), detP = -detP; j = find(p==i); p([i,j]) = p([j,i]); end
  end
else                                                 % symmetric B = I+sW*K*sW
  sW = sqrt(W); L = chol(eye(n)+sW*sW'.*K);             % Cholesky factor of B
  solveKiW = @(r) bsxfun(@times,L\(L'\bsxfun(@times,r,sW)),sW);
end
  
alpha = solveKiW(y_train-meanfunc(x_train, meanhyp));

% Estimate kernel output for covariance between training and test data
Ks = covfunc(x_train, x_test);
ms = meanfunc(x_test, meanhyp);                          % evaluate mean vector
mu = ms + Ks'*alpha(:,:);                               % conditional mean fs|f
