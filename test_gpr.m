function test_gpr

load /Volumes/UltraMax/BrainAGE_core/s8rp1_4mm_IXI547_r1840.mat
%load /Volumes/UltraMax/BrainAGE_core/s4rp1_4mm_UKB_1_r1700.mat

trend_degree = 1;
trend_method = 1;

[tmp, age_ind] = sort(age);
ind_test = age_ind(1:4:end);
ind_train = age_ind; ind_train(1:4:end) = [];

ind_roi = (1:size(Y,2)) > 0;

Y_train = Y(ind_train,ind_roi);
Y_test  = Y(ind_test,ind_roi);

age_train = age(ind_train);
age_test  = age(ind_test);

%----- ensure range 0..1
mn = min([Y_train(:); Y_test(:)]);
mx = max([Y_train(:); Y_test(:)]);
Y_train  = ((Y_train-mn)/(mx-mn));
Y_test   = ((Y_test-mn)/(mx-mn));

hyperparam = struct('mean', 100, 'lik', -1);

tic
[PredictedAge,vr] = cg_GPR3(Y_train, age_train, Y_test);
toc

BrainAGE = PredictedAge-age_test;
MAE_before = mean(abs(BrainAGE))
BrainAGE2 = apply_trend_correction(BrainAGE,age_test,trend_degree,trend_method);
MAE_after = mean(abs(BrainAGE2))

figure(11)
scatter(age_test,BrainAGE)

figure(12)
scatter(age_test,BrainAGE2)

tic
PredictedAge = cg_GPR(Y_train, age_train, Y_test);
toc


BrainAGE = PredictedAge-age_test;
MAE_before = mean(abs(BrainAGE))
BrainAGE2 = apply_trend_correction(BrainAGE,age_test,trend_degree,trend_method);
MAE_after = mean(abs(BrainAGE2))

figure(13)
scatter(age_test,BrainAGE)
figure(14)
scatter(age_test,BrainAGE2)

    
%-------------------------------------------------------------------------------
function [mu,vr] = cg_GPR3(x_train, y_train, x_test)

n_subj = numel(y_train);

xs = size(x_train,1);
if xs ~= n_subj
  x_train = x_train';
  x_test = x_test';
end

% Define covariance function
covfunc = @(x, y) x*y';

% Define prior mean function
meanfunc = @(x, mn) mn*ones(size(x,1),1);

% Set hyperparameters
meanhyp = 100;
likhyp  = log10(0.1); % noise variance

% Train GP
K = covfunc(x_train, x_train);

sn2 = exp(2*likhyp);                               % noise variance of likGauss
sW  = sqrt(1/sn2);
L   = chol(K/sn2 + eye(size(K,1)));

alpha = bsxfun(@times,L\(L'\bsxfun(@times,y_train-meanfunc(x_train, meanhyp),sW)),sW);

Ks = covfunc(x_train, x_test);

ms = meanfunc(x_test, meanhyp);                          % evaluate mean vector
mu = ms + Ks'*alpha(:,:);        

K_post = K + sn2*eye(size(K));
vr = diag(covfunc(x_test, x_test) - Ks'*((K_post\Ks)));

%-------------------------------------------------------------------------------
function mu = cg_GPR2(x_train, y_train, x_test)
%-------------------------------------------------------------------------------

n_subj = numel(y_train);

xs = size(x_train,1);
if xs ~= n_subj
  x_train = x_train';
  x_test = x_test';
end

% Define covariance function
covfunc = @(x, y) x*y';

% Define prior mean function
meanfunc = @(x, mn) mn*ones(size(x,1),1);

% Set hyperparameters
meanhyp = 100;
likhyp  = log10(0.1); % noise variance

% Train GP
K = covfunc(x_train, x_train);

sn2 = exp(2*likhyp);                               % noise variance of likGauss
W = ones(size(K,1),1)/sn2;

sW = sqrt(W); 
L = chol(eye(numel(W))+sW*sW'.*K);             % Cholesky factor of B
alpha = bsxfun(@times,L\(L'\bsxfun(@times,y_train-meanfunc(x_train, meanhyp),sW)),sW);

Ks = covfunc(x_train, x_test);

ms = meanfunc(x_test, meanhyp);                          % evaluate mean vector
mu = ms + Ks'*alpha(:,:);                               % conditional mean fs|f

%-------------------------------------------------------------------------------
function mu = cg_GPR_dropout(x_train, y_train, x_test, p_dropout)
%-------------------------------------------------------------------------------

n_subj = numel(y_train);

xs = size(x_train,1);
if xs ~= n_subj
  x_train = x_train';
  x_test = x_test';
end

n_voxel = size(x_train,2);

% Define covariance function with dropout
% Here we use simple linear covariance function
covfunc = @(x, y) x*y';

% Define prior (constant) mean function
meanfunc = @(x, mn) mn*ones(size(x,1),1);

% Set hyperparameters
meanhyp = 100;
likhyp  = log10(0.1); % noise variance

% number of Monte Carlo samples
if p_dropout
  n_mc = 10;
else
  n_mc = 1;
end

% Estimate kernel output for training data
K = covfunc(x_train, x_train);  
  
sn2 = exp(2*likhyp);                                % noise variance of likGauss
W   = ones(size(K,1),1)/sn2;

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
    % exclude random voxels w.r.t. dropout probability
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

%-------------------------------------------------------------------------------
function BrainAGE0 = apply_trend_correction(BrainAGE,age,trend_degree,trend_method)
%-------------------------------------------------------------------------------
% BrainAGE     - adjusted BrainAGE
% PredictedAge - adjusted PredictedAge
% Adjustment   - array of adjustments for indexed data

PredictedAge = BrainAGE + age;
if trend_degree < 0
  return
end

cc = corrcoef(BrainAGE,age);

% use BrainAGE for obtaining trend
if trend_method == 1
  Y = BrainAGE;
else
  % use predicted age for obtaining trend
  Y = PredictedAge;
end

% polynomial trend
if trend_degree > 0
  G = [ones(length(age),1) cg_polynomial(age,trend_degree)];
  G_indexed = [ones(length(age),1) cg_polynomial(age,trend_degree)];
% use just offset
else
  G = ones(length(age),1);
  G_indexed = ones(length(age),1);
end
  
% estimate beta only for indexed data (e.g. control subjects)
Beta = pinv(G_indexed)*Y;

% and remove effects for all data
GBeta = G*Beta;

% use BrainAGE for obtaining trend
if trend_method == 1
  PredictedAge  = PredictedAge  - GBeta;
  BrainAGE0 = PredictedAge - age;
else
  % use predicted age for obtaining trend
  PredictedAge  = (PredictedAge - Beta(1));
  if trend_degree == 1
    PredictedAge = PredictedAge/Beta(2);
  end
  BrainAGE0 = PredictedAge - age;
end
BrainAGE0 = BrainAGE0 - mean(BrainAGE0);
cc0 = corrcoef(BrainAGE0,age);
fprintf('Correlation to age before/after correction: %3.3f\t%3.3f\n',cc(1,2),cc0(1,2));