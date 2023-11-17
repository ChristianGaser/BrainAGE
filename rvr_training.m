function [dat,a] =  rvr_training(a,dat)
% modified spider function from @relvm_r/training.m

if(a.algorithm.verbosity>0)
  disp(['training ' get_name(a) '.... '])
end

t = get_y(dat);
if(size(t,2)>1)
    a = multi_reg(a);
    [dat,a] = train(a,dat);
    return;
end

if strcmp(a.child.ker,'poly')
  x = get_x(dat);
  % exclude NaNs
  x(isnan(x)) = 0;
  PHI = x*x' + 1;
  if a.child.poly > 1
    PHI = PHI.^a.child.poly;
  end
else
  fprintf('Use original spider code because kernel is not poly\n');
  PHI = calc(a.child,dat);
end

[N,M] = size(PHI);

w = zeros(M,1);
PHIt  = PHI'*t;

beta    = a.beta0;
alpha   = ones(M,1);
gamma   = ones(M,1);
nonZero = logical(ones(M,1));    
maxIts  = a.maxIts;

PRUNE_POINT = a.maxIts * (a.prune_point/100);
LAST_IT   = 0;
ALPHA_MAX   = 1e12;
MIN_DELTA_LOGALPHA  = 1e-2;    

for i=1:maxIts
    % 
    % Prune large values of alpha
    % 
%    nonZero = (alpha<ALPHA_MAX);
%?????????? is that better working?
    nonZero = (alpha<ALPHA_MAX & alpha >0);
    alpha_nz  = alpha(nonZero);
    w(~nonZero) = 0;
    M   = sum(nonZero);
    % Work with non-pruned basis
    % 
    PHI_nz  = PHI(:,nonZero);
    flag = 1;
    count = 0;

    Hessian = (PHI_nz'*PHI_nz);
    while flag
      count = count + 1;
      [U,flag]    = chol(Hessian*beta + diag(alpha_nz));
      
      % sometimes Hessian is not positive definite if beta is too large
      % mainly for 1st iteration (e.g. for large data sets)
      if flag
        if i == 1
          beta = 0.1;
        else
          beta = beta/10;
        end
%        fprintf('Hessian is not positive definite. Decrease beta parameter to %g after %d iterations.\n', beta,i);
      else
        count = 0;
      end
      % only allow 3 failed trials
      if count > 3
        break
      end
    end
    if flag
      fprintf('Terminating after %d iterations.\n',i)
      
      % return output arguments
      weights = w(nonZero);
      used  = find(nonZero);

      a.Xsv=get(dat,used);
      a.used=used;
      a.marginal=marginal;
      % !!! Note this difference!
      a.alpha=weights;
      a.beta=beta;
      a.gamma=gamma;
      a.trained=1;
      return
    else
      Ui    = inv(U);
      w(nonZero)  = (Ui * (Ui' * PHIt(nonZero)))*beta;
    end
    
    ED    = sum((t-PHI_nz*w(nonZero)).^2); % Data error
    betaED  = beta*ED;
    logBeta = N*log(beta);
    % Quick ways to get determinant and diagonal of posterior weight
    % covariance matrix 'SIGMA'
    logdetH = -2*sum(log(diag(Ui)));
    diagSig = sum(Ui.^2,2);
    % well-determinedness parameters
    gamma   = 1 - alpha_nz.*diagSig;
        
    % Compute marginal likelihood (approximation for classification case)
    marginal  = -0.5* (logdetH - sum(log(alpha_nz)) - ...
        logBeta + betaED + (w(nonZero).^2)'*alpha_nz);
    
    if(marginal==-Inf)
        break;
    end
    
    % Output info if requested and appropriate monitoring iteration
    if(a.algorithm.verbosity>1)
      if i>1
        fprintf('%5d> L = %.3f\t Gamma = %.2f (nz = %d)\t beta=%.3f\tmaxDAlpha=%g\n',...
            i, marginal, sum(gamma), sum(nonZero), beta, maxDAlpha);
      else
        fprintf('%5d> L = %.3f\t Gamma = %.2f (nz = %d)\t beta=%.3f\n',...
            i, marginal, sum(gamma), sum(nonZero), beta);
      end
    end
    if ~LAST_IT
        % 
        % alpha and beta re-estimation on all but last iteration
        % (we just update the posterior statistics the last time around)
        % 
        logAlpha    = log(alpha(nonZero));
        if i<PRUNE_POINT
            % MacKay-style update given in original NIPS paper
            alpha(nonZero)  = gamma ./ w(nonZero).^2;
        else
            % Hybrid update based on NIPS theory paper and AISTATS submission
            alpha(nonZero)  = gamma ./ (w(nonZero).^2./gamma - diagSig);
            alpha(alpha<=0) = inf; % This will be pruned later
        end
        anz   = alpha(nonZero);
        maxDAlpha = max(abs(logAlpha(anz~=0)-log(anz(anz~=0))));
        
        % Terminate if the largest alpha change is judged too small
        if maxDAlpha<MIN_DELTA_LOGALPHA
            LAST_IT = 1;
            if(a.algorithm.verbosity>0)
                fprintf('Terminating: max log(alpha) change is %g (<%g).\n', ...
                    maxDAlpha, MIN_DELTA_LOGALPHA);
            end
        end
        %
        % Beta re-estimate in regression (unless fixed)
        % 
        beta  = (N - sum(gamma))/ED;
    else
        % Its the last iteration due to termination, leave outer loop
        break;  % that's all folks!
    end
end
% Tidy up return values
weights = w(nonZero);
used  = find(nonZero);

if a.algorithm.verbosity>1
    fprintf('*\nHyperparameter estimation complete\n');
    fprintf('non-zero parameters:\t%d\n', length(weights));
    fprintf('log10-alpha min/max:\t%.2f/%.2f\n', ...
        log10([min(alpha_nz) max(alpha_nz)]));
end

a.Xsv=get(dat,used);
a.used=used;
a.marginal=marginal;
% !!! Note this difference!
a.alpha=weights;
a.beta=beta;
a.gamma=gamma;
a.trained=1;

if (~a.algorithm.do_not_evaluate_training_error)
  dat=test(a,dat);    
end
