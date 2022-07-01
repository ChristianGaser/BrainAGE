function [BrainAGE,EstimatedAge] = cg_BrainAGE(D)
%
% D.Y_test          - volume data for estimation
% D.age_test        - age of each volume 
% D.training_sample - training sample (e.g. {'IXI547','OASIS316'}
%     Healthy adults:
%         IXI547    - IXI-database
%         IXI410    - IXI-database (1:4:547 removed)
%         OASIS316  - OASIS
%         Ship1000  - SHIP
%     ADNI231Normal - ADNI Normal-sc-1.5T 
%Long_Australia2611 - PATH study Australia
%
%     Children data:
%         NIH394    - NIH objective 1
%         NIH876    - NIH release 4.0
%         NIH755    - NIH release 5.0
%
%     Children + adults:
%         fCONN772 - fcon-1000 (8-85 years)
%
% D.kernel          - kernel for rlvm_r approach (default kernel('poly',1))
% D.seg             - segmentation
%                     {'rp1'} use GM
%                     {'rp2'} use WM
%                     {'rp1,'rp2'} use both GM+WM
% D.hemi            - hemisphere for BrainAGE estimation
%                     'lh' only use left hemisphere
%                     'rh' only use right hemisphere
%                     'li' calculate lateralization index LI=(lh-rh)/(lh+rh)
%                     use both hemispheres as default
% D.res             - spatial resolution of data
% D.smooth          - smoothing size
% D.relnumber       - VBM release (e.g. '_r432')
% D.age_range       - age range of training data
%                     [0 Inf] use all data
%                     [50 80] use age range of 50..80
%                     if not defined use min/max of age of test data
% D.nuisance        - additionally define nuisance parameter for covarying out (e.g. gender)
% D.ind_train       - define indices of subjects used for training (e.g. limit the training to male subjects only)
% D.add_train       - add parameters to training data (e.g. sex)
% D.add_test        - add parameters to test data (e.g. sex)
% D.PCA             - apply PCA as feature reduction (default=1)
% D.PCA_method      - method for PCA
%                    'eig' Eigenvalue Decomposition of the covariance matrix (faster but less accurate, for compatibiliy)
%                    'svd' Singular Value Decomposition of X (the default)
% D.add_sex         - add sex to mapped data (default=0), only valid for D.PCA=1
% D.dir             - directory for databases and code
% D.spider_dir      - directory of spider
% D.verbose         - verbose level (default=1), set to "0" to suppress long outputs
% D.threshold_std   - all data with a standard deviation > D.threshold_std of mean covariance are excluded
%                     (after covarying out effects of age)
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_BrainAGE.m 2015-08-28 08:57:46Z gaser $


if ~isfield(D,'Y_test')
  error('D.Y_test not defined');
end

if ~isfield(D,'age_test')
  error('D.age_test not defined');
end

if ~isfield(D,'training_sample')
  error('D.training_sample not defined');
end

if ~isfield(D,'kernel')
  D.kernel = kernel('poly',1);
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

if ~isfield(D,'add_train')
  D.add_train = [];
  D.add_test  = [];
elseif ~isfield(D,'add_test')
  D.add_train = [];
  D.add_test  = [];
end

if ~isfield(D,'dir')
  D.dir = '/Volumes/UltraMax/BrainAGE_core';
end

if ~isfield(D,'threshold_std')
  D.threshold_std = Inf;
end

if ~isfield(D,'verbose')
  D.verbose = 1;
end

% add spider paths
if ~isfield(D,'spider_dir')
  D.spider_dir = '/Volumes/UltraMax/spider';
end

if ~exist('spider')
  addpath(D.spider_dir)
  spider_path;
end
eval(['addpath ' D.dir]);

if isfield(D,'training_sample')
  n_training_samples = length(D.training_sample);
else
  n_training_samples = 0;
end

Y_train    = [];
Y_train2   = [];
age_train  = [];
male_train = [];

% don't load training sample if test sample is the same (e.g. for n-fold validation)
if numel(D.train_array) == 1 && strcmp(D.train_array{1},D.data)
  load_training_sample = false;
  Y_train    = D.Y_test;
  age_train  = D.age_test;
  age        = D.age_test;
  male       = D.male_test;
  male_train = D.male_test;
else
  load_training_sample = true;
end

% load training sample(s)
for i = 1:n_training_samples
  if load_training_sample
    name = fullfile(D.dir, [D.smooth D.seg{1} '_' D.res 'mm_' D.training_sample{i} D.relnumber '.mat']);
    
    % check whether release number is already included in name
    if ~exist(name,'file')
      if contains(D.training_sample{i},'_r')
        name = fullfile(D.dir, [D.smooth D.seg{1} '_' D.res 'mm_' D.training_sample{i} '.mat']);
      end
    end
    
    if D.verbose > 1, fprintf('cg_BrainAGE: load %s\n',name); end
    load(name);
    Y_train    = [Y_train; single(Y)]; clear Y
    age_train  = [age_train; age];

    if ~exist('male','var')
      male = ones(size(age));
    end
    male_train = [male_train; male];
    
    if isfield(D,'hemi') && ~isempty(D.hemi)
      hemi_name = fullfile(D.dir, ['ind_' D.res 'mm.mat']);
      load(hemi_name);
    end

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
      Y_train2   = [Y_train2; single(Y)]; clear Y
    end

  end  
end
  
if isfield(D,'hemi')
  if strcmp(D.hemi,'lh')
    fprintf('Only use left hemisphere for BrainAGE estimation.\n');
    Y_train  = Y_train(:,ind_left);
    D.Y_test = D.Y_test(:,ind_left);
    if length(D.seg) > 1
      Y_train2 = Y_train2(:,ind_left);    
    end
  end
  if strcmp(D.hemi,'rh')
    fprintf('Only use right hemisphere for BrainAGE estimation.\n');
    Y_train = Y_train(:,ind_right);
    D.Y_test = D.Y_test(:,ind_right);
    if length(D.seg) > 1
      Y_train2 = Y_train2(:,ind_right);    
    end
  end
end

if length(D.seg) > 1 && load_training_sample
  Y_train = [Y_train Y_train2]; clear Y_train2
end

% apply comcat harmonization while preserving age effects
if  isfield(D,'comcat') && load_training_sample
  if length(D.comcat) ~= length(age_train)
    error('Size of site definition in D.comcat (n=%d) differs from sample size (n=%d)\n',length(D.comcat),length(age_train));
  end
  Y_train = cat_stat_comcat(Y_train, D.comcat, age_train, 0, 3, 1);
end

% use only indicated subjects (e.g. for gender-wise training)
if isfield(D,'ind_train')
  Y_train    = Y_train(D.ind_train,:);
  age_train  = age_train(D.ind_train);
  male_train = male_train(D.ind_train);
end

ind_age    = find(age_train >= D.age_range(1) & age_train <= D.age_range(2));
age_train  = age_train(ind_age);
male_train = male_train(ind_age);
Y_train    = Y_train(ind_age,:);

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
  if D.verbose, fprintf('%d subjects removed because their mean covariance was deviating more than %g standard deviation from mean.\n',length(ind_removed),D.threshold_std); end

  age_train(ind_removed) = [];
  male_train(ind_removed) = [];
  Y_train(ind_removed,:) = [];
end

D.male_train = male_train;

% give warning if age range differs by more than two years
if (min(age_train)-min(D.age_test) > 2) || (max(D.age_test)>max(age_train) > 2)
  fprintf('Warning: Defined age range of training sample (%g..%g years) does not fit to real age range of sample %g..%g\n',...
      min(age_train),max(age_train),min(D.age_test),max(D.age_test));
end

if D.verbose
  fprintf('%d subjects used for training (age %3.1f..%3.1f years)\n',length(age_train),min(age_train),max(age_train));
  fprintf('Mean age\t%g (SD %g) years\nMales/Females\t%d/%d\n',mean(age_train),std(age_train),sum(male_train),length(age_train)-sum(male_train));
end

%--- KERNEL -> RVR
s = relvm_r(D.kernel);

% get rid of excessive output
s.algorithm.verbosity = D.verbose;
s.maxIts = 250;

% use smaller beta for larger training samples to increase stability
if isfield(s,'beta0') && (numel(age_train) > 1000)
  s.beta0 = 0.1;
end

%----- ensure range 0..1
mn = min([Y_train(:); D.Y_test(:)]);
mx = max([Y_train(:); D.Y_test(:)]);
Y_train  = ((Y_train-mn)/(mx-mn));
D.Y_test = ((D.Y_test-mn)/(mx-mn));

% size of training and test sample
n_train = size(Y_train, 1);
n_test  = size(D.Y_test, 1);

%---PCA---
if D.PCA
  % if defined limit number of PCA components
  if D.PCA > 1
    n_PCA = D.PCA;
  else % otherwise use n-1 components (minimum of voxel or subjects)
    n_PCA = min(size(Y_train)) - 1;
  end

  [mapped_train, mapping] = cg_pca(Y_train, n_PCA,D.PCA_method);

  clear Y_train
  mapped_test = (D.Y_test - repmat(mapping.mean, [size(D.Y_test, 1) 1])) * mapping.M;
  clear mapping
  
  mapped_train = double(mapped_train);
  mapped_test  = double(mapped_test);
  
  % ensure range 0..1
  mn = min([mapped_train(:);mapped_test(:)]);
  mx = max([mapped_train(:);mapped_test(:)]);
  mapped_train = (mapped_train - mn)/(mx - mn);
  mapped_test  = (mapped_test - mn)/(mx - mn);
else
  mapped_train = double(Y_train); clear Y_train
  mapped_test  = double(D.Y_test);
end

% add information about sex to mapped data
if isfield(D,'add_sex') && D.add_sex > 0
  if ~isempty(D.male_train) && ~isempty(D.male_test)
    mapped_train = [mapped_train D.male_train];
    mapped_test  = [mapped_test D.male_test];
  else
    fprintf('Information about sex not found in mat-file.\n');
  end
end

% add additional parameters to mapped data
if ~isempty(D.add_train) && ~isempty(D.add_test)
  if size(D.add_train, 1) ~= n_train
    error('Data size between training data (n=%d) and additional parameter (n=%d) differs.\n');
  end
  if size(D.add_test, 1) ~= n_test
    error('Data size between test data (n=%d) and additional parameter (n=%d) differs.\n');
  end
  mapped_train = [mapped_train D.add_train];
  mapped_test  = [mapped_test D.add_test];
end

d  = data(mapped_train, age_train);
t1 = data(mapped_test, D.age_test);

%---RVR---
% I use a modified rvr_training because the original @relvm_r/training
% function has issues with stability and I have added an additional break
try
%    evalc('[est,model] = rvr_training(s,d);');
  [est,model] = rvr_training(s,d);
  
  % Test
  if 1
  pred1 = test(model,t1);
  EstimatedAge = pred1.X;
  else
    EstimatedAge = rvr_test(model,t1);
  end
  
  % sometimes chol function in rvr_training fails and we obtain zero values
  % for EstimatedAge that will be replaced by NaNs
  if std(EstimatedAge) == 0
    BrainAGE = NaN(size(D.age_test));
    EstimatedAge = NaN(size(D.age_test));
  else
    if isfield(D,'define_cov')
      BrainAGE = EstimatedAge;
    else
      BrainAGE = EstimatedAge-D.age_test;
    end
  end
  
catch
  BrainAGE = NaN(size(D.age_test));
  EstimatedAge = NaN(size(D.age_test));
  fprintf('\n-----------------------------------------\n');
  fprintf('ERROR: Training failed.\n');
  fprintf('-----------------------------------------\n\n');
  return
end

if ~isempty(D.nuisance)
  G = [ones(length(D.age_test),1) D.nuisance];
  if D.verbose, fprintf('Remove effect of nuisance parameter\n'); end
 
  % estimate beta
  Beta = pinv(G)*BrainAGE;

  % and remove effects
  BrainAGE = BrainAGE - G*Beta;
  EstimatedAge = EstimatedAge - G*Beta;
end

return
