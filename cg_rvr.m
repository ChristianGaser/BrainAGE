function [varargout] = cg_rvr(varargin)
% Optimisation for Relevance Vector Regression
%
% [w,alpha,beta,ll] = cg_rvr(Phi,t)
% Phi   - MxM matrix derived from kernel function of vector pairs
% t     - the values to be matched
% w     - weights
% alpha - 1/variance for the prior part of the model
% beta  - 1/variance for the likelihood part of the model
% ll    - the negative log-likelihood.
%
% [w,alpha,beta,nu,ll]=cg_rvr(K,t,opt)
% K     - a cell-array of MxM dot-product matrices.
% t     - the values to be matched
% opt   - either 'Linear' or 'Gaussian RBF'
%         'Linear'       is for linear regression models, where
%                        the optimal kernel is generated by
%                        [nu(1)*K{1} + nu(1)*K{2}... ones(size(K{1},1),1)]
%         'Gaussian RBF' is for regression using Gaussian radial basis
%                        functions.  The kernel is generated from
%                        P1  = nu(1)*K{1} + nu(1)*K{2} ... ;
%                        P2  = repmat(diag(P1) ,1,size(P1,2)) +...
%                              repmat(diag(P1)',size(P1,1),1) - 2*P1;
%                        Phi = exp([-0.5*P2 ones(size(P1,1),1)]);
% w     - weights
% alpha - 1/variance for the prior part of the model
% beta  - 1/variance for the likelihood part of the model
% nu    - parameters that convert the dot-product matrices into
%         a kernel matrix (Phi).
% ll    - the negative log-likelihood.
%
% The first way of calling the routine simply optimises the
% weights.  This involves estimating a restricted maximum
% likelihood (REML) solution, which maximises P(alpha,beta|t,Phi).
% Note that REML is also known as Type II Maximum Likelihood
% (ML-II). The ML-II solution tends towards infinite weights for
% some the regularisation terms (i.e. 1/alpha(i) approaches 0).
% The appropriate columns are removed from the model when
% this happens.
%
% The second way of calling the routine also estimates additional
% input scale parameters as described in Appendix C of Tipping (2001).
% This method is much slower, as a full optimisation for the scale
% parameters is done after each update of the alphas and beta.
%
% see: http://research.microsoft.com/mlp/RVM/relevance.htm
%
% Refs:
% The Relevance Vector Machine.
% In S. A. Solla, T. K. Leen, and K.-R. Müller (Eds.),
% Advances in Neural Information Processing Systems 12,
% pp.  652-658. Cambridge, Mass: MIT Press.
%
% Michael E. Tipping
% Sparse Bayesian Learning and the Relevance Vector Machine
% Journal of Machine Learning Research 1 (2001) 211-244
%________________________________________________________________________
% Copyright (C) 2011 Wellcome Department of Imaging Neuroscience & Machine Learning & Neuroimaging Laboratory

% Written by John Ashburner
% $Id: cg_rvr.m 542 2012-05-20 10:48:51Z cphillip $

[varargout{1:nargout}]=regression0(varargin{:});
return;
%__________________________________________________________________________

%__________________________________________________________________________
function [w,alpha,beta,ll]=regression0(Phi,t)
[N,M]  = size(Phi);
if N==M
    Phi = [Phi ones(N,1)];
elseif M~=N+1,
    error('Phi must be N x (N+1)');
end;
scale             = sqrt(sum(sum(Phi(1:N,1:N).^2))/N^2);
scale             = [ones(N,1)*scale ; 1];
Phi               = Phi/spdiags(scale,0,numel(scale),numel(scale));
alpha             = ones(size(Phi,2),1)/N;
%beta             = N/sum((t-mean(t)).^2);
beta              = 1e6;
[w,alpha,beta,ll] = rvr1a(Phi,t,alpha,beta);
alpha             = [alpha(1)*ones(N,1) ; alpha(2)];
[w,alpha,beta,ll] = rvr2a(Phi,t,alpha,beta);
w                 = w./scale;
alpha             = alpha.*scale.^2;
return;
%__________________________________________________________________________

%__________________________________________________________________________
function [w,alpha,beta,ll] = rvr1a(Phi,t,alpha,beta)
% This function is not actually used
%spm_chi2_plot('Init','ML-II (non-sparse)','-Log-likelihood','Iteration');
[N,M]   = size(Phi);
ll      = Inf;
PP      = Phi'*Phi;
Pt      = Phi'*t;
for subit=1:10,
    alpha_old = alpha;
    beta_old  = beta;

    % E-step
    S         = inv(PP*beta + spdiags([ones(N,1)*alpha(1) ; alpha(2)],0,N+1,N+1));
    w         = S*(Pt*beta);

   % figure(3); plot(t,Phi*w,'.'); drawnow;

    tmp       = t-Phi*w;
    ll        = ...
         -0.5*log(alpha(1))*N-0.5*log(alpha(2))-0.5*N*log(beta)-0.5*logdet(S)...
         +0.5*tmp'*tmp*beta + 0.5*sum(w.^2.*[repmat(alpha(1),N,1) ; alpha(2)])...
         +0.5*(M-N)*log(2*pi);
  %  if subit>1, spm_chi2_plot('Set',ll); end;
    %fprintf('%g\n',ll);

    % M-step
    ds        = diag(S);
    dfa1      = sum(ds(1:N))*alpha(1);
    dfa2      = sum(ds(N+1))*alpha(2);
    alpha(1)  = max(N-dfa1,eps)/(sum(w(1:N).^2)   +eps);
    alpha(2)  = max(1-dfa2,eps)/(sum(w(N+1).^2)   +eps);
    beta      = max(dfa1+dfa2-1,eps)/(sum((Phi*w-t).^2)+eps);

    % Convergence
    if max(max(abs(log((alpha+eps)./(alpha_old+eps)))),log(beta/beta_old)) < 1e-9,
        break;
    end;
end;
%spm_chi2_plot('Clear');
return;
%__________________________________________________________________________

%__________________________________________________________________________
function [w,alpha,beta,ll]=rvr2a(Phi,t,alpha,beta)
%spm_chi2_plot('Init','ML-II (sparse)','-Log-likelihood','Iteration');
[N,M] = size(Phi);
nz    = true(M,1);

PP    = Phi'*Phi;
Pt    = Phi'*t;

for subit=1:200,
    th         = min(alpha)*1e9;
    nz         = alpha<th;
    alpha(~nz) = th*1e9;
    alpha_old  = alpha;
    beta_old   = beta;

    % E-step
    S         = inv(PP(nz,nz)*beta + diag(alpha(nz)));
    w         = S*Pt(nz)*beta;

  %  figure(3); plot(t,Phi(:,nz)*w,'.'); drawnow;

    tmp = t-Phi(:,nz)*w;
    ll  = ...
        -0.5*sum(log(alpha(nz)+1e-32))-0.5*N*log(beta+1e-32)-0.5*logdet(S)...
        +0.5*tmp'*tmp*beta + 0.5*sum(w.^2.*alpha(nz))...
        +0.5*(sum(nz)-N)*log(2*pi);
%    if subit>0, spm_chi2_plot('Set',ll); end;
    %fprintf('%g\t%g\n',ll,exp(mean(log(alpha)))/beta);

    % M-step
    gam       = 1 - alpha(nz).*diag(S);
    alpha(nz) = max(gam,eps)./(w.^2+1e-32);
    beta      = max(N-sum(gam),eps)./(sum((Phi(:,nz)*w-t).^2)+1e-32);

    % Convergence
    if max(max(abs(log((alpha(nz)+eps)./(alpha_old(nz)+eps)))),log(beta/beta_old)) < 1e-6*N,
        break;
    end;
end;
w(nz)  = w;
w(~nz) = 0;
w      = w(:);
%spm_chi2_plot('Clear');
%__________________________________________________________________________


%__________________________________________________________________________
function [ld,C] = logdet(A)
A  = (A+A')/2;
C  = chol(A);
d  = max(diag(C),eps);
ld = sum(2*log(d));

