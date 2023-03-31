function test_rvr

addpath('/Volumes/UltraMax/spider')
spider_path;

trend_degree = 2;

%--- KERNEL -> RVR
s = relvm_r(kernel('poly',1));

% get rid of excessive output
s.algorithm.verbosity = 1;
s.beta0 = 100;

load /Volumes/UltraMax/BrainAGE_core/s4rp1_8mm_IXI547_r1840.mat

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

d  = data(Y_train, age_train);
t1 = data(Y_test, age_test);

n = [1:0.1:2];
n = 5;
p = n./2.6;

x = zeros(numel(n),numel(p));

tic
for in=1:length(n)
   for ip=1:length(p)
    
    sn = n(in);
    sp = p(ip);

  [m_post,m_pred,c_pred] = gpr(Y_train,age_train,Y_test,sn,sp);
  
  PredictedAge = m_pred;

  BrainAGE2 = PredictedAge-age_test;
  fprintf('%3.3f\t%3.2f\t%3.2f\n',mean(abs(BrainAGE2)),sn,sp)
  x(in,ip) = mean(abs(BrainAGE2));
end
end
toc


%imagesc(x)
%colorbar
[est,model] = rvr_training(s,d);
pred1 = test(model,t1);

EstimatedAge = pred1.X;
BrainAGE1 = EstimatedAge-age_test;
mean(abs(BrainAGE1))

addpath('gpml-v4.2')
gp_startup

meanfunc = @meanConst;
covfunc = @covRQiso;
likfunc = @likGauss;

hyperparam = struct('mean', 5, 'cov', [3 3 0.1], 'lik', -1);
covfunc = @covLINiso;
hyperparam = struct('mean', 100, 'cov', 0, 'lik', -1);

if 0
R = 1;
[N,D] = size(Y_train);
covfunc = { 'covADD',{1:R,'covSEiso'} };  
hyperparam.cov = [ log(ones(1,2*D)), log(ones(1,R))];    
end

tic
[PredictedAge, Var_PredictedAge] = gp(hyperparam, @infExact, meanfunc, covfunc, likfunc, Y_train, age_train, Y_test);
toc
BrainAGE3 = PredictedAge-age_test;
MAE = mean(abs(BrainAGE3))

tic
hyp2 = minimize(hyperparam, @gp, -100, @infExact, meanfunc, covfunc, likfunc, Y_train, age_train)
toc

[PredictedAge, Var_PredictedAge] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, Y_train, age_train, Y_test);
BrainAGE3 = PredictedAge-age_test;
MAE = mean(abs(BrainAGE3))


figure(11)
scatter(PredictedAge,Var_PredictedAge)
figure(12)
scatter(age_test,PredictedAge-age_test)




function [mu, s2] = gp_regression(x_train, y_train, x_test, hyp)

% Define covariance function
covfunc = @(x, y, hyp) hyp(1)^2 * exp(-0.5*hyp(2)^2*(pdist2(x',y')/hyp(3)^2).^2);

% Define likelihood function
likfunc = @(f, y, hyp) sum(log(normpdf(y, f, sqrt(exp(hyp)))));

% Define prior mean function
meanfunc = @(x) zeros(size(x));

% Set hyperparameters
meanhyp = [];
likhyp = log(0.1); % noise variance
covhyp = [1, 1, 1]; % hyperparameters for the covariance function

% Train GP
K = covfunc(x_train, x_train, covhyp);
L = chol(K + exp(likhyp)*eye(size(K)));

alpha = solve_chol(L, y_train-meanfunc(x_train), true);
beta = solve_chol(L, ones(length(x_train), 1), true);
f = meanfunc(x_test) + covfunc(x_train, x_test, covhyp).' * alpha;
v = covfunc(x_test, x_test, covhyp) - covfunc(x_train, x_test, covhyp).' * solve_chol(L, covfunc(x_train, x_test, covhyp), false);
s2 = diag(v);
mu = f;


function [x] = solve_chol(L, b, transpose_flag)
% Solves the system L'*L*x=b or L*L'*x=b, depending on the transpose_flag.
if transpose_flag
    x = L.' \ (L \ b);
else
    x = L \ (L.' \ b);
end

function K = apx(hyp,cov,x,opt)

function [ldB2,solveKiW,dW,dldB2,L,triB] = ldB2_exact(W,K,dK)
  isWneg = any(W<0); n = numel(W);
  if isWneg                  % switch between Cholesky and LU decomposition mode
    A = eye(n) + bsxfun(@times,K,W');                     % asymmetric A = I+K*W
    [L,U,P] = lu(A); u = diag(U);         % compute LU decomposition, A = P'*L*U
    signU = prod(sign(u));                                           % sign of U
    detP = 1;               % compute sign (and det) of the permutation matrix P
    p = P*(1:n)';
    for i=1:n                                                     % swap entries
      if i~=p(i), detP = -detP; j = find(p==i); p([i,j]) = p([j,i]); end
    end
    if signU~=detP     % log becomes complex for negative values, encoded by inf
      ldB2 = Inf;
    else          % det(L) = 1 and U triangular => det(A) = det(P)*prod(diag(U))
      ldB2 = sum(log(abs(u)))/2;
    end                                            % compute inverse if required
    if nargout>1, Q = U\(L\P); solveKiW = @(r) bsxfun(@times,W,Q*r); end
    if nargout>4, L = -diag(W)*Q; end                              % overwrite L
  else                                                 % symmetric B = I+sW*K*sW
    sW = sqrt(W); L = chol(eye(n)+sW*sW'.*K);             % Cholesky factor of B
    ldB2 = sum(log(diag(L)));                                    % log(det(B))/2
    solveKiW = @(r) bsxfun(@times,solve_chol(L,bsxfun(@times,r,sW)),sW);
    if nargout>2, Q = bsxfun(@times,1./sW,solve_chol(L,diag(sW))); end
  end
  if nargout>2
    dW = sum(Q.*K,2)/2;            % d log(det(B))/2 / d W = diag(inv(inv(K)+W))
    triB = trace(Q);                                      % triB = trace(inv(B))
    dldB2 = @(varargin) ldB2_deriv_exact(W,dK,Q, varargin{:});     % derivatives
  end
  
  
function [post nlZ dnlZ] = infGaussLik(hyp, mean, cov, lik, x, y)

% Exact inference for a GP with Gaussian likelihood.
%
% Compute a parametrization of the posterior, the negative log marginal
% likelihood and its derivatives w.r.t. the hyperparameters. The function takes
% a specified covariance function (see covFunctions.m) and likelihood function
% (see likFunctions.m), and is designed to be used with gp.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2018-08-22.
%                                      File automatically generated using noweb.
%
% See also INFMETHODS.M, APX.M.

[n, D] = size(x);
[m,dm] = feval(mean{:}, hyp.mean, x);           % evaluate mean vector and deriv
sn2 = exp(2*hyp.lik); W = ones(n,1)/sn2;            % noise variance of likGauss

[K,dK] = feval(cov{:},hyp.cov,x);     % covariance matrix and dir derivative
[ldB2,solveKiW,dW,dhyp,post.L] = ldB2_exact(W,K,dK); % obtain functionality depending on W

alpha = solveKiW(y-m);
post.alpha = alpha;                            % return the posterior parameters
post.sW = sqrt(W);                              % sqrt of noise precision vector
if nargout>1                               % do we want the marginal likelihood?
  nlZ = (y-m)'*alpha/2 + ldB2 + n*log(2*pi*sn2)/2;    % -log marginal likelihood
  if nargout>2                                         % do we want derivatives?
    dnlZ = dhyp(alpha); dnlZ.mean = -dm(alpha);
    dnlZ.lik = -sn2*(alpha'*alpha) - 2*sum(dW)/sn2 + n;
  end
end


function [m_post,m_pred,c_pred] = gpr(Y,x,yt,sn,sp)

% GZ 22/3/2012
% 
% note: here we apply the 'inverted Rasmussen notation' to allow multiple 
% applications in neuroimaging given brains Y and corresponding labels x 
% (columns from covariate matrix X)
%
% goal: makes predictions for unseen brains yt, given a multivariate input dataset 
% Y and corresponding labels x, uses the dual implementation to calculate
% the mean and the variance of the predictive distribution
%      
%    p(x_est|yt,Y,x)
% 
% Makes the following simplifying assumtions: 

% Firstly, the feature space embedding the identity, i.e. from the rasmussen 
% perspective phi(y)=y and thus the feature space is identical to the n 
% dimensional inputspace.
%  
% Secondly, we assume a gaussian prior with Sp=sp*speye(n) covariance matrix for
% the weightings which simplifies further things a lot. This could be
% basically replaced by a spatial prior.
% 
% Thirdly, the model of y is a linear function of the inputs i.e.
%
%   f(x)=x'w + e, e ~ N(0,sn^2)
% 
% The sigma^2 must be provided. 

% INPUT
% Y'  - ell x n brain datamat (is internally transposed) with n dim
%       inputspace
%       observations in rows
% x   - ell x 1 covariate vector (e.g. age) used as dependent variable for
%       regression 
% yt  - testbrain for prediction
% sigmasq - error variance of the univariate model error 

%disp 'linear bayesian regression with gaussian prior';
%clc;
% training data size
n_subj = numel(x);

[~, ys] = size(Y);
if ys ~= n_subj
  Y = Y';
  yt = yt';
end

% prior covariance
% kernel matrix i.e. K=phi(Y)'*SIGMA_p*phi(Y)=Y'*diag(sp*ones(n,1))*Y
K = Y'*(sp^2*Y);

% dual variables
alpha = (sp^2*Y)/(K+sn^2*eye(n_subj));

% posterior distribution p(w|Y,x) ~ N(sn^(-2)*inv(A)*Y*x,inv(A))
% using dual expression: sn^(-2)*inv(A)*Y=Sp*Y*inv(K+sn^2*eye(ell))=alpha
% mean
m_post = alpha*x;

% evaluation of the predictive distribution
% mean
m_pred = yt'*m_post;

% covariance 
c_pred=diag(yt'*(sp^2*yt)-yt'*alpha*Y'*(sp^2*yt));


