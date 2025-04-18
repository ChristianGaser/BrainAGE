function [BrainAGE, PredictedAge, D] = BA_gpr(D)
% [BrainAGE, PredictedAge, D] = BA_gpr_ui(D)
% BrainAGE estimation using Gaussian Process Regression (GPR)
%
% D.Y_test          - data for estimation
% D.age_test        - age of each data set 
%
% D.hyperparam      - GPR hyperparameters (.mean and .lik)
% D.seg             - segmentation
%                     {'rp1'} use GM
%                     {'rp2'} use WM
%                     {'rp1,'rp2'} use both GM+WM
% D.res             - spatial resolution of data
% D.smooth          - smoothing size
% D.relnumber       - CAT12 release (e.g. '_r1840')
% D.age_range       - age range of training data
%                     [0 Inf] use all data
%                     [50 80] use age range of 50..80
%                     if not defined use min/max of age of test data
% D.comcat          - If data are acquired at different sites (e.g. using different scanners or sequences) we can
%                     harmonize data using ComCAT. A vector with coding of the scanners is required for the test data (EXPERIMENTAL!).
% D.ind_train       - define indices of subjects used for training (e.g. limit the training to male subjects only)
% D.RVR             - use old RVR method (Spider toolbox necessary, where the toolbox path can be defined using D.spider_dir)
% D.PCA             - apply PCA as feature reduction (default=1)
% D.PCA_method      - method for PCA
%                    'eig' Eigenvalue Decomposition of the covariance matrix (faster but less accurate, for compatibiliy)
%                    'svd' Singular Value Decomposition of X (the default)
% D.dir             - directory for databases and code
% D.p_dropout       - Dropout probability to randomly exclude voxels/data points to implement an uncertainty-aware approach using a 
%                     Monte-Carlo Dropout during inference. That means that during testing, voxels are randomly dropped out according 
%                     to the dropout probabilities. This process is repeated multiple times, and each time, the model produces 
%                     a different output. By averaging these outputs, we can obtain a more robust prediction and estimate the model's 
%                     uncertainty in its predictions. A meaningful dropout probability is 0.1, which means that 10% of the data points 
%                     are excluded. The default is 0.
% D.verbose         - verbose level (default=1), set to "0" to suppress long outputs
% D.threshold_std   - all data with a standard deviation > D.threshold_std of mean covariance are excluded
%                     (after covarying out effects of age)
% D.eqdist          - options for age and sex equalization between test and train
% D.eqdist.weight   - vector of size 2 that allows to weight the cost function for age and sex equalization
% D.eqdist.range    - matrix 2 x 2 which defines the age range and sex range for equalization
% D.eqdist.tol      - vector of size 2 that defines tolerance between mean value of
%                     age_test and age_train and male_test and male_train
% D.nuisance        - additionally define nuisance parameter for covarying out (e.g. gender)
% D.mask            - define mask to additionally estimate regional BrainAGE values 
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if ~isfield(D,'Y_test')
  error('D.Y_test not defined');
end

if ~isfield(D,'age_test')
  error('D.age_test not defined');
end

if ~isfield(D,'training_sample')
  error('D.training_sample not defined');
end

if ~isfield(D,'hyperparam')
  D.hyperparam = struct('mean', 100, 'lik', -1);
end

if ~isfield(D,'RVR')
  D.RVR = 0;
end

if ~isfield(D,'PCA')
  D.PCA = 1;
end

if ~isfield(D,'PCA_method')
  D.PCA_method = 'svd';
end

if ~iscell(D.training_sample)
  D.training_sample = cellstr(D.training_sample);
end

if ~isfield(D,'seg')
  D.seg = 'rp1';
end

if ~iscell(D.seg)
  D.seg = cellstr(D.seg);
end

if ~isfield(D,'res')
  D.res = 8;
end

if ~isfield(D,'smooth')
  D.smooth = 's8';
end

if ~isfield(D,'relnumber')
  D.relnumber = '_r432';
end

if ~isfield(D,'age_range')
  D.age_range = [min(D.age_test) max(D.age_test)];
end

if ~isfield(D,'nuisance')
  D.nuisance = [];
end

if ~isfield(D,'p_dropout')
  D.p_dropout = 0;
end

if ~isfield(D,'dir')
  D.dir = '/Volumes/UltraMax/BrainAGE';
end

if ~isfield(D,'threshold_std')
  D.threshold_std = Inf;
end

if ~isfield(D,'verbose')
  D.verbose = 1;
end

if isfield(D,'training_sample')
  n_training_samples = length(D.training_sample);
else
  n_training_samples = 0;
end

Y_train    = [];
Y_train2   = [];
age_train  = [];
male_train = [];

if D.RVR
  if ~isfield(D,'spider_dir')
    D.spider_dir = '/Volumes/UltraMax/spider';
  end

  addpath(D.spider_dir)
  spider_path;
  
  % KERNEL -> RVR
  s = relvm_r(kernel('poly',1));
  
  % get rid of excessive output
  s.algorithm.verbosity = D.verbose;
  s.maxIts = 250;

  % use smaller beta for larger training samples to increase stability
  if isfield(s,'beta0') && (numel(age_train) > 1000)
    s.beta0 = 0.1;
  end
end

Y_test = D.Y_test;

% don't load training sample if test sample is the same (e.g. for k-fold validation)
if numel(D.train_array) == 1 && strcmp(D.train_array{1},D.data)
  load_training_sample = false;
  Y_train    = Y_test;
  age_train  = D.age_test;
  age        = D.age_test;
  
  if isfield(D,'male_test') && ~isempty(D.male_test)
    male = D.male_test;
    male_train = D.male_test; 
  else
    male = ones(size(age));
    male_train = ones(size(age)); 
  end  
else
  load_training_sample = true;  
end

% load training sample(s)
for i = 1:n_training_samples
  if load_training_sample
    name0 = [D.smooth D.seg{1} '_' D.res 'mm_' D.training_sample{i} D.relnumber '.mat'];
    name = fullfile(D.dir, name0);
    
    % check whether release number is already included in name
    if ~exist(name,'file')
      if contains(D.training_sample{i},'_r')
        name = fullfile(D.dir, [D.smooth D.seg{1} '_' D.res 'mm_' D.training_sample{i} '.mat']);
      end
    end
    
    if D.verbose > 1, fprintf('BA_gpr_core: load %s\n',name); end

    % search first if the file exists in current folder
    if ~exist(name) && exist(name0)
      load(name0)
    else load(name); end
    
    Y_train    = [Y_train; single(Y)]; clear Y
    age_train  = [age_train; age];
        
    if ~exist('male','var') || (exist('male','var') && isempty(male))
      male = ones(size(age));
    end
    male_train = [male_train; male];
    
    % consider additional segmentations
    if length(D.seg) > 1
      for j = 2:length(D.seg)
        name = fullfile(D.dir, [D.smooth D.seg{j} '_' D.res 'mm_' D.training_sample{i} D.relnumber '.mat']);
        
        % check whether release number is already included in name
        if ~exist(name,'file')
          if contains(D.training_sample{i},'_r')
          name = fullfile(D.dir, [D.smooth D.seg{j} '_' D.res 'mm_' D.training_sample{i} '.mat']);
          end
        end
        load(name);
      end
      Y_train2 = [Y_train2; single(Y)]; clear Y
    end

  end  
end

if isfield(D,'mask')
  [m,n] = size(Y_train);
  if n ~= numel(D.mask)
    error('Size mismatch for mask. There are %d entries necessary',n);
  end
  if ~islogical(D.mask)
    D.mask = D.mask > 0.5;
  end
  
  % restrict data for training and test to mask area
  Y_train  = Y_train(:,D.mask);
  Y_test = Y_test(:,D.mask);
end

if length(D.seg) > 1 && load_training_sample
  Y_train = [Y_train Y_train2]; clear Y_train2
end

% apply comcat harmonization while preserving age effects
if  isfield(D,'comcat')
  if length(D.comcat) ~= length(D.age_test)
    error('Size of site definition in D.comcat (n=%d) differs from test sample size (n=%d)\n',length(D.comcat),length(D.age_test));
  end
  fprintf('Apply ComCat for %d sites\n',numel(unique(D.comcat)));
  Y_test = cat_stat_comcat(Y_test, D.comcat,[],D.age_test, 0, 1, 0);
end

% use only indicated subjects (e.g. for gender-wise training or k-fold cross-validation)
% skip that option for minimizing hyperparameters because we need the full training sample
if isfield(D,'ind_train')
  Y_train    = Y_train(D.ind_train,:);
  age_train  = age_train(D.ind_train);
  try
    male_train = male_train(D.ind_train);
  catch
    male_train = nan(size(D.ind_train));
  end
end

ind_age    = find(age_train >= D.age_range(1) & age_train <= D.age_range(2));
if isfield(D,'k_fold') && D.k_fold > 0 && numel(ind_age) < numel(age_train)
  fprintf('Warning: Limiting age range is not allowed for k-fold validation.\n');
else
  age_train  = age_train(ind_age);
  if ~isempty(male_train), male_train = male_train(ind_age); end
  Y_train    = Y_train(ind_age,:);
end

% remove subjects where standard deviation of mean covariance > D.threshold_std
if ~isinf(D.threshold_std)
  Y = Y_train;
  n_subjects = size(Y,1);
  
  % remove age effects from Y to reliably estimate covariance without its influence
  Ymean = repmat(mean(Y), [n_subjects 1]);
  G = age_train - mean(age_train);
  Y = Y - G*(pinv(G)*Y) + Ymean; clear Ymean

  % calculate covariance matrix
  YpY = (Y*Y')/n_subjects; clear Y

  % normalize YpY
  d      = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
  dd     = d*d';
  YpY    = YpY./(dd+eps);
  t      = find(abs(YpY) > 1); 
  YpY(t) = YpY(t)./abs(YpY(t));
  YpY(1:n_subjects+1:end) = sign(diag(YpY));

  % extract mean correlation for each data set
  mean_cov = zeros(n_subjects,1);
  for i=1:n_subjects
    % extract row for each subject
    cov0 = YpY(i,:);

    % remove cov with its own
    cov0(i) = [];
    mean_cov(i) = mean(cov0);
  end

  % sort files
  [mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
  threshold_cov = mean(mean_cov) - D.threshold_std*std(mean_cov);
  n_thresholded = min(find(mean_cov_sorted < threshold_cov));

  ind_removed = find(mean_cov < threshold_cov);
  if D.verbose > 1, fprintf('%d subjects removed because their mean covariance was deviating more than %g standard deviation from mean.\n',length(ind_removed),D.threshold_std); end

  age_train(ind_removed)  = [];
  male_train(ind_removed) = [];
  Y_train(ind_removed,:)  = [];
end

% if at least tolerance is defined for age and sex equalization we can select 
% train data to better match test and train distribution
if isfield(D,'eqdist') && isfield(D.eqdist,'tol') 
  
  n_before = length(age_train);
  
  % check whether we already have to estimate assignement or we have
  % multiple training sample
  if ~isfield(D.eqdist,'sel') || numel(D.train_array) > 1
    if numel(D.train_array) == 1 && strcmp(D.train_array{1},D.data)
      fprintf('Use of eqdist flag is not allowed for k-fold validation\n');
      return
    end

    % if not defined use age range from given D.age_range and don't limit sex parameter
    if ~isfield(D.eqdist,'range'), D.eqdist.range = [D.age_range; 0 1]'; end

    % create sample for age/sex for test as reference and train as source
    if ~isfield(D,'male_test'), D.male_test = round(rand(size(D.age_test))); end
    ind_groups = [];
    for o = 1:numel(D.ind_groups)
      ind_groups = [ind_groups; D.ind_groups{o}];
    end
    sample_ref = [D.age_test(ind_groups) D.male_test(ind_groups)];
    sample_src = [age_train male_train];

    % find assignement using Hungarian method
    fprintf('Estimate assignement to match distribution of train data to test data\n');
    [sample_sel, sel] = BA_equalize_distribution(sample_ref, sample_src, D.eqdist);
    fprintf('Selected %d of %d data sets for training after age/sex equalization.\n',numel(sel),numel(age_train));

    % save selection to skip time consuming assignement for next runs
    D.eqdist.sel = sel;
  else
    sel = D.eqdist.sel;
  end
  % get age, sex, and Ytrain back
  age_train  = age_train(sel);
  male_train = male_train(sel);
  Y_train = Y_train(sel,:);
  
  if D.verbose && (n_before-length(age_train)) > 0, fprintf('%d subjects removed to equalize distribution of age and sex.\n',n_before-length(age_train)); end
  
  D.male_train = male_train;
end


% give warning if age range differs by more than two years
if (D.age_range(1)-min(D.age_test)) > 2 || (max(D.age_test)-D.age_range(2)) > 2
  fprintf('\n********************************************************************\n'); 
  warning('Warning: Defined age range of training sample (%3.1f..%3.1f years) differs from real age range of sample %g..%g\n',...
      min(age_train),max(age_train),min(D.age_test),max(D.age_test));
  fprintf('********************************************************************\n'); 
end

% save output
D.age_train = age_train;
D.male_train = male_train;

% estimate range using training data only and scale to 0..1
mn = min(Y_train(:));
mx = max(Y_train(:));
Y_train  = (Y_train-mn)/(mx-mn);
Y_test = (Y_test-mn)/(mx-mn);

% PCA using training data only
if D.PCA 
  % if defined limit number of PCA components
  if D.PCA > 1
    n_PCA = D.PCA;
  else % otherwise use n-1 components (minimum of voxel or subjects)
    n_PCA = min(size(Y_train)) - 1;
  end

  [mapped_train, mapping] = BA_pca(Y_train,n_PCA,D.PCA_method);

  clear Y_train
  mapped_test = (Y_test - repmat(mapping.mean, [size(Y_test, 1) 1])) * mapping.M;
  if ~D.p_dropout, clear mapping; end
  
  % we have to use double format for GPR
  mapped_train = double(mapped_train);
  mapped_test  = double(mapped_test);

  % ensure range 0..1
  mn = min(mapped_train(:));
  mx = max(mapped_train(:));
  mapped_train = (mapped_train - mn)/(mx - mn);
  mapped_test  = (mapped_test - mn)/(mx - mn);
  
  mapping.mn = mn;
  mapping.mx = mx;
  
else
  % we have to use double format for GPR
  mapped_train = double(Y_train); clear Y_train
  mapped_test  = double(Y_test);
end

% Regression using GPR
if ~D.RVR
  if D.PCA && D.p_dropout
    PredictedAge = BA_gpr_core(mapped_train, age_train, mapped_test, ...
            D.hyperparam.mean, D.hyperparam.lik, D.p_dropout, mapping, D.Y_test);
  else
    PredictedAge = BA_gpr_core(mapped_train, age_train, mapped_test, ...
            D.hyperparam.mean, D.hyperparam.lik);
  end
  BrainAGE = PredictedAge-D.age_test;
else
  % Regression using RVR

  d  = data(mapped_train, age_train);
  t1 = data(mapped_test, D.age_test);
  
  try
    [~,model] = rvr_training(s,d);
  
    pred1 = test(model,t1);
    PredictedAge = pred1.X;
  catch
    fprintf('RVR-prediction failed.\n');
    PredictedAge = nan(size(D.age_test));
  end

  % sometimes chol function in rvr_training fails and we obtain zero values
  % for EstimatedAge that will be replaced by NaNs
  if std(PredictedAge) == 0
    BrainAGE = NaN(size(D.age_test));
    PredictedAge = NaN(size(D.age_test));
  else
    BrainAGE = PredictedAge-D.age_test;
  end
end


if ~isempty(D.nuisance)
  G = [ones(length(D.age_test),1) D.nuisance];
  if D.verbose > 1, fprintf('Remove effect of nuisance parameter\n'); end
 
  % estimate beta
  Beta = pinv(G)*BrainAGE;

  % and remove effects
  BrainAGE = BrainAGE - G*Beta;
  PredictedAge = PredictedAge - G*Beta;
end

