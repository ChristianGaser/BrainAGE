function [BrainAGE, BrainAGE_unsorted, BrainAGE_all, D, age] = cg_BrainAGE_GPR_ui(D)
% [BrainAGE, BrainAGE_unsorted,BrainAGE_all,D] = cg_BrainAGE_GPR_ui(D)
% user interface for BrainAGE estimation
%
% D.train_array     - cell array of training samples
%     Healthy adults:
%         IXI547     - IXI-database
%         IXI410     - IXI-database (1:4:547 removed)
%         OASIS316   - OASIS
%         OASIS_3381 - OASIS3
%         Ship1000   - SHIP (internal use only)
%     ADNI231Normal  - ADNI Normal-sc-1.5T 
%       UKB_x_r1700  - UKB-data sample x with release r1700 
%
%     Children data:
%         NIH394     - NIH objective 1
%         NIH876     - NIH release 4.0
%         NIH755     - NIH release 5.0
%
%     Children + adults:
%         fCONN772  - fcon-1000 (8-85 years)
%
% D.hyperparam      - GP hyperparameters
% D.data            - test sample for BrainAGE estimation
% D.seg_array       - segmentation
%                     {'rp1'} use GM
%                     {'rp2'} use WM
%                     {'rp1,'rp2'} use GM or WM
%                     {'rp1+rp2'} use both GM+WM
% D.res_array       - spatial resolution of data as char or cell (can be a cell array to try different values: {'4','8'})
% D.smooth_array    - smoothing size as char or cell (can be a cell array to try different values: {'s4','s8'})
% D.hemi            - hemisphere for BrainAGE estimation
%                     'lh' only use left hemisphere
%                     'rh' only use right hemisphere
%                     'li' calculate lateralization index LI=(lh-rh)/(lh+rh)
%                     use both hemispheres as default
% D.relnumber       - VBM release (e.g. '_r432')
% D.age_range       - age range of training data
%                     [0 Inf] use all data (default)
%                     [50 80] use age range of 50..80
%                     if not defined use min/max of age of test data
% D.nuisance        - additionally define nuisance parameter for covarying out (e.g. gender)
% D.ind_adjust      - define indices for adjusting data according to trend defined with D.trend_degree
%                     usually this is the control group that is used in order to adjust the data
% D.site_adjust     - If data are acquired at different sites (e.g. using different scanners or sequences) the trend 
%                     correction should be estimated and applied for each site seperately. A vector with coding of the scanners is required.
%                     If this parameter is empty this ensures that no site-specific adjustment is made even for multiple
%                     training data that are combined with a "+". 
% D.comcat          - If data are acquired at different sites (e.g. using different scanners or sequences) we can
%                     harmonize data using ComCAT. A vector with coding of the scanners is required (EXPERIMENTAL!).
% D.ind_groups      - define indices of groups, if D.ind_adjust is not given, then the first group of this index will
%                     be used for adjusting data according to trend defined with D.trend_degree
% D.ind_train       - define indices of subjects used for training (e.g. limit the training to male subjects only)
% D.trend_degree    - estimate trend with defined order using healthy controls and apply it to all data (set to -1 for skipping trend correction)
% D.trend_method    - use different methods for estimating trend:
%                     1 use BrainAGE for trend correction (default)
%                     2 use predicted age for trend correction (as used in Cole et al. 2018)
% D.PCA             - apply PCA as feature reduction (default=1), values > 1 define number of PCA components
% D.PCA_method      - method for PCA
%                    'eig' Eigenvalue Decomposition of the covariance matrix (faster but less accurate, for compatibiliy)
%                    'svd' Singular Value Decomposition of X (the default)
% D.k_fold          - k-fold validation if training and test sample are the same or only one is defined (10-fold as default)
%                     Common approach for k-fold is to divide the sample into k parts and to use the
%                     larger part (n-n/k) for training and the remaining part (n/k) for testing.
%                     For huge data sets such as UKB, the training data are getting too large and we 
%                     switch the selection: we now use the smaller part (n/k) for training and the
%                     larger part (n-n/k) for testing) by simply defining a negative value for k_fold
%                     The additional entries (a multiple of k-1) that result from that approach (by factor k-1) are finally averaged.
% D.k_fold_TPs      - definition of time points for k-fold validation to ensure that multiple time points of one subject are not mixed 
%                     between test and training data (only necessary to define for longitudinal data and k-fold validation)
% D.k_fold_reps     - Number of repeated k-fold cross-validation
% D.k_fold_rand     - As default the age values for the training sample is sorted and every k-th data is selected for training to minimize age 
%                     differences between training and test data. With k_fold_rand you can set the seed for the random number generator.
% D.ensemble        - ensemble method to combine different models
%                     0 - Majority voting: use model with lowest MAE
%                     1 - Weighted GLM average: use GLM estimation to estimate model weights to minimize MAE
%                     2 - Median: use median to weight different models (as default)
%                     3 - GLM: use GLM estimation to maximize variance to a group or a regression parameter (EXPERIMENTAL!)
%                     4 - Stacking: use GPR to combine models (EXPERIMENTAL!, only works with k_fold validation)
%                     5 - Weighted Average: (average models with weighting w.r.t. MAE)
% D.contrast        - define contrast to maximize group differences (use only if D.ensemble=3) (e.g. [1 -1])
%                     D.contrast can be also a vector which is used to maximize variance between BrainAGE and this parameter.
% D.dir             - directory for databases and code
% D.verbose         - verbose level (default=1), set to "0" to suppress long outputs
% D.threshold_std   - all data with a standard deviation > D.threshold_std of mean covariance are excluded
%                     (after covarying out effects of age)
% D.eqdist          - options for age and sex equalization between test and train
% D.eqdist.weight   - vector of size 2 that allows to weight the cost function for age and sex equalization
% D.eqdist.range    - matrix 2 x 2 which defines the age range and sex range for equalization
% D.eqdist.tol      - vector of size 2 that defines tolerance between mean value of
%                     age_test and age_train and male_test and male_train
% D.eqdist.debug    - print debug info if set (default 0)
% D.corr            - additionally define parameter that can be correlated to BrainAGE if only one group is given
% D.define_cov      - optionally define continous parameter that should be used instead of age for more general use 
%                     of GPR not only limited to BrainAGE
% D.style           - plot-style: 1: old style with vertical violin-plot; 2: new style with horizontal density plot
% D.groupcolor      - matrix with (group)-bar-color(s), use jet(numel(data)) or other color functions (nejm by default)
%
% Parameter search
% ---------------
% Some selected parameters can be also defined as ranges to try different parameter settings.
% Examples:
% D.trend_degree = 2;
% D.threshold_std = [Inf];
% D.age_range = [20 50];
% D.res_array    = {'4','8'};   
% D.smooth_array = {'s4','s8'}; 
% D.train_array = {'IXI547','OASIS316+ADNI231Normal'}; 
% D.ensemble = 1; % minimize MAE
% D.data = 'Your_Sample';
%
% Output
% -------
% BrainAGE            - BrainAGE values sorted by group definitions in D.ind_groups 
% BrainAGE_unsorted   - unsorted (originally ordered) BrainAGE values
% BrainAGE_all        - array of BrainAGE values for all models sorted by group definitions in D.ind_groups
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_BrainAGE_GPR_ui.m 4 2015-08-28 08:57:46Z gaser $

global min_hyperparam

addpath(fullfile(fileparts(mfilename('fullpath')),'gpml-v4.2'))
gp_startup

if ~isfield(D,'trend_degree')
  D.trend_degree = 2;
end

if ~isfield(D,'trend_method')
  D.trend_method = 1;
end

if ~isfield(D,'age_range')
  D.age_range = [0 Inf];
end

if ~isfield(D,'hyperparam')
  D.hyperparam = struct('mean', 100, 'lik', -1);
end

% Set the mean function, covariance function and likelihood function
% Take meanConst, covRQiso and likGauss as default
if ~isfield(D,'meanfunc')
  D.meanfunc = @meanConst;
end

if ~isfield(D,'covfunc')
  D.covfunc = @covLIN;
end

if ~isfield(D,'likfunc')
  D.likfunc = @likGauss;
end

if ~isfield(D,'minimize_hyperparam')
  D.minimize_hyperparam = false;
end

if ~isfield(D,'style')
  style = 2;
else
  style = D.style;
end

if ~isfield(D,'threshold_std')
  D.threshold_std = Inf;
end

if ~isfield(D,'PCA')
  D.PCA = 1;
end

if ~isfield(D,'PCA_method')
  D.PCA_method = 'svd';
end

if D.trend_method > 1 && D.trend_degree > 1
  D.trend_degree = 1;
  fprintf('Only use linear trend correction for method that uses predicted age for obtaining trend correction.\n');
end

% this is just for compatbility with older scripts
if isfield(D,'seg') && ~isfield(D,'seg_array')
  D.seg_array = D.seg;
  D = rmfield(D,'seg');
end

% this is just for compatbility with older scripts
if isfield(D,'training_sample') && ~isfield(D,'train_array')
  if numel(D.training_sample) > 1
    error('Please use new syntax and use D.train_array instead.');
  end
  D.train_array = D.training_sample;
  D = rmfield(D,'training_sample');
end

% convert to cell if necessary
if ~iscell(D.seg_array)
  D.seg_array = cellstr(D.seg_array);
end

% array with different smoothing sizes
if ~iscell(D.smooth_array)
  D.smooth_array = cellstr(D.smooth_array);
end

% array with different spatial resolutions
if ~iscell(D.res_array)
  D.res_array = cellstr(D.res_array);
end

% verbose level
if ~isfield(D,'verbose')
  D.verbose = 1;
end

% fill the missing field if neccessary
if ~isfield(D,'data')
  D.data = D.train_array{1};
end

% set default for k_fold_reps 
if ~isfield(D,'k_fold_reps')
  D.k_fold_reps = 1;
end

if isfield(D,'k_fold_rand') && D.k_fold_reps > 1
  error('D.k_fold_rand cannot be used together with D.k_fold_reps because repeated k-fold would always use the same random numbers without variations.');
end

if iscell(D.data)
  D.data = char(D.data);
end

% consider old syntax and name
if isfield(D,'n_fold') && ~isfield(D,'k_fold')
  D.k_fold = D.n_fold;
end

ensemble_str = {'Majority Voting (model with lowest MAE)',...
                'Weighted GLM Average (GLM for weighting models to minimize MAE)',...
                'Median of all models',...
                'GLM for weighting models to maximize variance w.r.t. contrast vector',...
                'Stacking: GPR for combining models',...
                'Weighted Average: (average models with weighting w.r.t. MAE)'};

if isfield(D,'define_cov')
  if ~isempty(strfind(D.data,'+'))
    error('D.define_cov cannot be used for multiple training data.');
  end
  if isfield(D,'n_data')
    if numel(D.define_cov) ~= D.n_data
      error('D.define_cov has different size (%d) than data (%d).',numel(D.define_cov),D.n_data);
    end
  else
    % assume that D.define_cov is correctly defined and we can obtain
    % data size form that
    D.n_data = numel(D.define_cov);
    age = D.define_cov;
  end
  if D.trend_degree > -1
    D.trend_degree = -1;
    fprintf('Disable trend correction because this cannot be used for other non-BrainAGE parameters\n');
  end
end

% load first data set to get data size
if ~isfield(D,'n_data')
  ind_plus = strfind(D.seg_array{1},'+');
  if ~isempty(ind_plus)
    seg_array = D.seg_array{1}(1:ind_plus-1);
  else
    seg_array = D.seg_array{1};
  end

  % find potential "+" indicating to combine training sample
  ind_plus = strfind(D.data,'+');
  if ~isempty(ind_plus)
    age0 = [];
    ind_plus = [0 ind_plus length(D.data)+1];
    n_train = numel(ind_plus)-1;
    for l=1:n_train
      load([D.smooth_array{1} seg_array '_' D.res_array{1} 'mm_' D.data(ind_plus(l)+1:ind_plus(l+1)-1) D.relnumber],'age');
      age0 = [age0; age]; 
    end
    age = age0;
    D.n_data = numel(age);
  else
    load([D.smooth_array{1} seg_array '_' D.res_array{1} 'mm_' D.data D.relnumber],'age');
    D.n_data = numel(age);
  end
end

if ~isfield(D,'ind_groups')
  D.ind_groups{1} = 1:D.n_data;
end

if ~isfield(D,'site_adjust')
  D.site_adjust = ones(D.n_data,1);
end

if ~isfield(D,'ind_adjust')
  D.ind_adjust = D.ind_groups{1};
end

if isfield(D,'run_kfold') && ~isfield(D,'train_array')
  D.train_array = {D.data};
end

if isfield(D,'weighting') && ~isfield(D,'ensemble')
  D.ensemble = D.weighting;
  fprintf('This option is deprecated. Use the option ''ensemble'' instead.\n');
  return
end

% check whether contrast was defined for ensemble=3
if isfield(D,'ensemble')
  if D.ensemble(1) == 3 && ~isfield(D,'contrast')
    error('D.contrast has to be defined.');
  end
end

% print some parameters
if ~isfield(D,'run_kfold')
  res = []; seg = []; smo = [];
  for i = 1:numel(D.res_array)
    res = [res D.res_array{i} ' '];
  end
  for i = 1:numel(D.smooth_array)
    smo = [smo D.smooth_array{i} ' '];
  end
  for i = 1:numel(D.seg_array)
    seg = [seg D.seg_array{i} ' '];
  end
  fprintf('--------------------------------------------------------------\n');
  fprintf('Data:         \t%s\nResolution:   \t%s\nSmoothing:    \t%s\nSegmentation:  \t%s\nThreshold-Std:\t%d\n',...
    [D.data D.relnumber],res,smo,seg,D.threshold_std);
  if isfield(D,'train_array')
    tra = [];
    for i = 1:numel(D.train_array)
      tra = [tra D.train_array{i} ' '];
    end
    fprintf('Training-Data:\t%s\n',tra);
  end
  if isfield(D,'ensemble')
    fprintf('Model-Weight: \t%d\n',D.ensemble);
  end
  if isfield(D,'k_fold')
    fprintf('k-Fold:       \t%d\n',D.k_fold);
  end
  fprintf('PCA:          \t%d (method: %s)\n',D.PCA,D.PCA_method);
  fprintf('Trend correct:\t%d\n',D.trend_degree);
  fprintf('Age-Range:    \t%g-%g\n',D.age_range(1),D.age_range(2));
  fprintf('--------------------------------------------------------------\n');
end

% run k-fold validation if no data field is given or validation with k_fold is defined
if ((~isfield(D,'data') || ~isfield(D,'train_array')) || isfield(D,'k_fold')) && ~isfield(D,'run_kfold')
  
  if isfield(D,'ensemble') && D.ensemble(1) == 3 && numel(D.ind_groups) < 1
    error('Weighting=3 cannot be used within k-fold validation with more than one group.');
  end

  ind_adjust = D.ind_adjust;
    
  % use 10-fold as default
  if ~isfield(D,'k_fold')
    D.k_fold = 10;
  end
    
  ind_all  = [];
  age_all  = [];
  BA_all   = [];
    
  if D.k_fold < 1
    D.k_fold = -D.k_fold;
    inverse_k_fold = 1;
  else
    inverse_k_fold = 0;
  end

  % ensure that this field is always defined and set to ones by default
  if ~isfield(D,'k_fold_TPs')
    D.k_fold_TPs = ones(D.n_data,1);
  end
  
  % number of time points
  n_TPs = max(D.k_fold_TPs);
  
  % for longitudinal data only
  if n_TPs > 1
    % find order of time point definition
    % offset_TPs = 1 -> alternating order (e.g. 1 2 1 2 1 2)
    % offset_TPs > 1 -> consecutive order (e.g. 1 1 1 2 2 2)
    offset_TPs = find(diff(D.k_fold_TPs), 1 );
    
    ns = [];
    for i=1:n_TPs
      ns = [ns sum(D.k_fold_TPs==i)];
    end
    
    if any(diff(ns))
      error('Time points should all have same size if you apply k-fold to longitudinal data.');
    end
    
    if offset_TPs == 1
      fprintf('Longitudinal data with %d time points with alternating order found.\n',n_TPs);
    else
      fprintf('Longitudinal data with %d time points with consecutive order found.\n',n_TPs);
    end

    % sort only data for TP1
    [~, ind_age] = sort(age(D.k_fold_TPs==1));
    
    if offset_TPs == 1
      ind_age = ind_age * n_TPs - (n_TPs-1);
    end
    
    for i=1:n_TPs-1
      ind_age = [ind_age; ind_age+i*offset_TPs];
    end
  else
    [~, ind_age] = sort(age);
  end

  min_hyperparam = cell(numel(D.res_array),numel(D.smooth_array),numel(D.seg_array),numel(D.train_array));

  for rep = 1:D.k_fold_reps
    
    % control random number generation to get always the same seeds
    if isfield(D,'k_fold_rand')
      rng(D.k_fold_rand)
      ind_age = randperm(numel(age))';
    end
    
    % use random indexing of age for repeated k-fold corss-validation
    if D.k_fold_reps > 1 && ~isfield(D,'k_fold_rand')
      if exist('rng','file') == 2
        rng('default')
        rng(rep)
      else
        rand('state',rep);
      end
      ind_age = randperm(numel(age))';
    end
  
    for j = 1:D.k_fold
      % indicate that validation is running and will not be called in nested loops
      D.run_kfold = j;
      D.run_repetition = rep;

      ind_fold0 = j:D.k_fold:D.n_data;

      % try to use similar age distribution between folds
      % by using sorted age index
      ind_test = ind_age(ind_fold0)';

      % build training sample using remaining subjects
      ind_train = ind_age';
      ind_train(ind_fold0) = [];

      % Common approach for k-fold is to divide the sample into k parts and to use the
      % larger part (n-n/k) for training and the remaining part (n/k) for testing.
      % For huge data sets such as UKB, the training data are getting too large and we 
      % switch the selection: we now use the smaller part (n/k) for training and the
      % larger part (n-n/k) for testing)
      if inverse_k_fold
        tmp_train = ind_train;
        ind_train = ind_test;
        ind_test  = tmp_train;
      end

      % I know this should never happen, but be absolutely sure we check
      % whether there is some overlap between training and test data
      n_overlaps = sum(ismember(ind_train,ind_test));
      if n_overlaps
        fprintf('WARNING: There is an overlap of %d subjects between training and test data.\n',n_overlaps)
      end

      % build indices for training and test
      ind_train_array{j} = ind_train;
      ind_test_array{j}  = ind_test;

      % collect age and indces in the order w.r.t. the folding
      age_all       = [age_all; age(ind_test)];
      ind_all       = [ind_all ind_test];

      % prepare cg_BrainAGE_GPR_ui parameters
      D.ind_groups  = {ind_test};
      D.ind_adjust  = ind_test;
      D.ind_train   = ind_train;

      % call nested loop
      [BA_fold_all, ~, ~, D] = cg_BrainAGE_GPR_ui(D);

      if j == 1 && rep == 1
        BA_all = zeros(D.n_data,D.k_fold,D.k_fold_reps,size(BA_fold_all,2));
      end

      % we keep entries for each loop because there might be some overlapping
      % entries if inverse_k_fold is used
      BA_all(ind_test,j,rep,:) = BA_fold_all;
    end
  end
  
  % we have to estimate the mean by using the sum and dividing by the
  % actual numbers of entries (~=0)
  n_entries = sum(BA_all(:,:,1)~=0,2);
  BA_all = squeeze(sum(BA_all,2))./n_entries;
  if any(n_entries>1)
    if min(n_entries) == max(n_entries)
      fprintf('There are %d overlapping entries that were averaged.\n',max(n_entries)); 
    else
      fprintf('There are %d-%d overlapping entries that were averaged.\n',max(n_entries),max(n_entries)); 
    end
  end

  % we use the mean over all repetitions
  if D.k_fold_reps > 1
    BA_all = squeeze(mean(BA_all,2));
  end
  
  D.ind_adjust = ind_adjust; % rescue original ind_adjust
  
  % go through different ensembles if defined
  D0 = D;
  
  D.MAE = [];
    
  if isfield(D,'ensemble')
    BA_unsorted_weighted = [];

    for m=1:numel(D.ensemble)
      D0.ensemble = D.ensemble(m);
      
      % apply initial trend correction to all models separately, but only for GPR
      % approach
      if D.trend_degree >= 0 && D0.ensemble == 4
        for i=1:size(BA_all,2)
          BA_all(:,i) = apply_trend_correction(BA_all(:,i),age,D,0);
        end
      end
      
      [~, PredictedAge_unsorted_weighted] = ensemble_models(BA_all,age,D0,ind_test_array,ind_train_array);

      BA_unsorted_weighted0 = PredictedAge_unsorted_weighted-age;
      
      % apply final trend correction to weighted model
      if D.trend_degree >= 0
        [BA_unsorted_weighted0,~,Adjustment] = apply_trend_correction(BA_unsorted_weighted0,age,D);
      end

      MAE_weighted = mean(abs(BA_unsorted_weighted0));
      D.MAE = [D.MAE MAE_weighted];
      if isfield(D,'define_cov')
        cc = corrcoef(BA_unsorted_weighted0,age);
      else
        cc = corrcoef(BA_unsorted_weighted0+age,age);
      end
      fprintf('\n===========================================================\n'); 
      fprintf(ensemble_str{D0.ensemble+1}); fprintf('\n');
      if ~isfield(D,'define_cov')
        fprintf('Overall weighted MAE for %d-fold = %g (ensemble=%d)\n',D.k_fold,MAE_weighted,D0.ensemble);
      end
      fprintf('Overall weighted correlation for %d-fold = %g\n',D.k_fold,cc(1,2));
      fprintf('============================================================\n\n'); 
      BA_unsorted_weighted = [BA_unsorted_weighted BA_unsorted_weighted0];
    end
  else
    BA_unsorted_weighted = BA_all;

    % apply trend correction
    if D.trend_degree >= 0
      [BA_unsorted_weighted,~,Adjustment] = apply_trend_correction(BA_unsorted_weighted,age,D);
    end
    
    % only print performance for single model
    if size(BA_unsorted_weighted,2) == 1
      D.MAE = mean(abs(BA_unsorted_weighted));
      if isfield(D,'define_cov')
        ind = ~isnan(age);
        cc = corrcoef(BA_unsorted_weighted(ind),age(ind));
      else
        cc = corrcoef(BA_unsorted_weighted+age,age);
      end
      fprintf('\n===========================================================\n'); 
      if ~isfield(D,'define_cov')
        fprintf('Overall MAE for %d-fold = %g\n',D.k_fold,D.MAE);
      end
      fprintf('Overall correlation for %d-fold = %g\n',D.k_fold,cc(1,2));
      fprintf('============================================================\n\n'); 
    end
  end
  
  BrainAGE_unsorted = BA_unsorted_weighted;
  BrainAGE = BrainAGE_unsorted;
  
  if nargout > 2
    % if site_adjust is empty we don't apply site adjustment
    if isempty(D.site_adjust)
      site_adjust = ones(D.n_data,1);
    else
      site_adjust = D.site_adjust;
    end
    
    % apply trend correction only for non-weighted data because otherwise
    % it's already done
    if D.trend_degree >= 0 && ~isfield(D,'ensemble')
      % apply trend correction for each site separately
      for i = 1:max(site_adjust)
        ind_site = find(site_adjust == i);
        BA_all(ind_site,:) = BA_all(ind_site,:) - Adjustment{i};
      end
    end      
    BrainAGE_all = BA_all;
  end
  
  return
end

% check whether additional fields for weighted BA are available
multiple_BA = numel(D.smooth_array) > 1 || numel(D.seg_array) > 1 || numel(D.res_array) > 1;

% prepare output
BA              = [];
BA_unsorted     = [];
PredictedAge_unsorted     = [];
EA              = [];

if ~exist('min_hyperparam','var')
  min_hyperparam = cell(numel(D.res_array),numel(D.smooth_array),numel(D.seg_array),numel(D.train_array));
end

orig_hyperparam = D.hyperparam;

% go through all resolutions, smoothing sizes, segmentations and training samples
for i = 1:numel(D.res_array)
  for j = 1:numel(D.smooth_array)
    for k = 1:numel(D.seg_array)
      for q = 1:numel(D.train_array)
        
        % select current data
        D.res = D.res_array{i};
        D.smooth = D.smooth_array{j};
        seg = D.seg_array{k};
        training_sample = D.train_array{q};

        % remove old D.training_sample field if exist
        if isfield(D,'training_sample')
          D = rmfield(D,'training_sample');
        end

        % remove old D.seg field if exist
        if isfield(D,'seg')
          D = rmfield(D,'seg');
        end
        
        % find potential "+" indicating to combine training sample
        ind_plus = strfind(training_sample,'+');
        if ~isempty(ind_plus)
          ind_plus = [0 ind_plus length(training_sample)+1];
          n_train = numel(ind_plus)-1;
          for l=1:n_train
            D.training_sample{l} = training_sample(ind_plus(l)+1:ind_plus(l+1)-1);
          end
        else
          D.training_sample{1} = training_sample;
        end

        % find potential "+" indicating to combine segmentations
        ind_plus = strfind(seg,'+');
        if ~isempty(ind_plus)
          ind_plus = [0 ind_plus length(seg)+1];
          n_seg = numel(ind_plus)-1;
          for l=1:n_seg
            D.seg{l} = seg(ind_plus(l)+1:ind_plus(l+1)-1);
          end
        else
          n_seg = 1;
          D.seg{1} = seg;
        end
        
        % clear old male parameter        
        if exist('male','var')
          clear male
        end
        
        % find potential "+" indicating to combine data
        ind_plus = strfind(D.data,'+');
        if ~isempty(ind_plus)
          Y0    = [];
          age0  = [];
          male0 = [];
          name0 = [];
          site_adjust = [];
          ind_plus = [0 ind_plus length(D.data)+1];
          n_train = numel(ind_plus)-1;
          for l=1:n_train
            name = [D.smooth D.seg{1} '_' D.res 'mm_' D.data(ind_plus(l)+1:ind_plus(l+1)-1) D.relnumber];
            if D.verbose > 1, fprintf('cg_BrainAGE_GPR_ui: load %s\n',name); end
            load(name);
            name0 = [name0 '+' name]; 
            age0  = [age0; age]; 
            male0 = [male0; male]; 
            Y0    = [Y0; single(Y)]; clear Y
            site_adjust = [site_adjust; l*ones(size(age))];  
          end
          age  = age0;
          male = male0;
          name = name0(2:end); % remove leading '+'
          Y    = Y0; clear Y0
          
          % create D.site_adjust if not already defined for more than one site
          if max(D.site_adjust) == 1
            D.site_adjust = site_adjust;
          elseif ~isempty(D.site_adjust)
            if i==1 && j==1 && k==1 && q==1
              fprintf('\n-----------------------------------------------------------------\n'); 
              fprintf('Please ensure that site-specific adjustment is correctly defined also for each training sample!\n');
              fprintf('\n-----------------------------------------------------------------\n'); 
            end
          end
        else
          name = [D.smooth D.seg{1} '_' D.res 'mm_' D.data D.relnumber];
          if D.verbose > 1, fprintf('cg_BrainAGE_GPR_ui: load %s\n',name); end
          load(name);
        end
        
        n_data = size(Y,1);
        if D.n_data ~= n_data
          fprintf('\n-----------------------------------------------------------------\n'); 
          fprintf('Data size differs for %s (%d vs. %d)',name,D.n_data,n_data);
          fprintf('\n-----------------------------------------------------------------\n'); 
          return
        end
    
        if D.verbose, fprintf('\n-----------------------------------------------------------------\n%s\n',name); end
        D.Y_test = single(Y); clear Y
        if isfield(D,'define_cov')
          age = D.define_cov;
          male = ones(size(age));
        end
        D.age_test = age;
        
        if exist('male','var')
          D.male_test = male;
        end
    
        if ~isfield(D,'age_range')
          D.age_range = [min(D.age_test) max(D.age_test)];
        end
        
        % use additional segmentation if defined
        for l = 2:n_seg
          % find potential "+" indicating to combine data
          ind_plus = strfind(D.data,'+');
          if ~isempty(ind_plus)
            Y0    = [];
            ind_plus = [0 ind_plus length(D.data)+1];
            n_train = numel(ind_plus)-1;
            for m=1:n_train
              name = [D.smooth D.seg{l} '_' D.res 'mm_' D.data(ind_plus(m)+1:ind_plus(m+1)-1) D.relnumber];
              if D.verbose > 1, fprintf('cg_BrainAGE_GPR_ui: load %s\n',name); end
              load(name);
              Y0    = [Y0; single(Y)]; clear Y
            end
            D.Y_test = [D.Y_test Y0]; clear Y0
          else
            name = [D.smooth D.seg{l} '_' D.res 'mm_' D.data D.relnumber];
            if D.verbose > 1, fprintf('cg_BrainAGE_GPR_ui: load %s\n',name); end
            load(name);
            D.Y_test = [D.Y_test single(Y)]; clear Y
          end
        end

        if D.verbose, fprintf('\n'); end
        
        % apply comcat harmonization while preserving age effects
        if isfield(D,'comcat')
          if length(D.comcat) ~= length(D.age_test)
            error('Size of site definition in D.comcat (n=%d) differs from sample size (n=%d)\n',...
              length(D.comcat),length(D.age_test));
          end
          D.Y_test = cat_stat_comcat(D.Y_test, D.comcat, [], D.age_test, 0, 3, 1);
        end
        
        if ~isfield(D,'ind_groups')
          D.ind_groups = {1:length(D.age_test)};
        end
        
        n_groups = numel(D.ind_groups);
        
        if ~isfield(D,'groupcolor')
          groupcolor = nejm;
        else
          groupcolor = D.groupcolor;
        end
    
        if ~isfield(D,'nuisance')
          D.nuisance = [];
        end
        
        % build index for test data
        ind_test = [];
        for o = 1:n_groups
          if size(D.ind_groups{o},1) < size(D.ind_groups{o},2)
            D.ind_groups{o} = D.ind_groups{o}';
          end
          ind_test = [ind_test; D.ind_groups{o}];
        end
                                

        if ~isfield(D,'age_range')
          D.age_range = [min(D.age_test(ind_test)) max(D.age_test(ind_test))];
        else
          if numel(D.age_range) ~=2
            error('Age range has to be defined by two values (min/max)');
          end
        end
                      
        % estimate lateralization index LI
        calc_LI = 0;
        if isfield(D,'hemi')
          if strcmp(D.hemi,'li')
            calc_LI = 1;
          end
        end
        
        % minimize hyperparameter and start minimization only for the first
        % run, otherwise use minimized hyperparameter
        if D.minimize_hyperparam
          if (isfield(D,'run_kfold') && D.run_kfold == 1)
            D0 = D;
            D.ind_train = 1:numel(D.age_test);
            D.hyperparam = orig_hyperparam;
            min_hyperparam{i,j,k,q}  = cg_BrainAGE_GPR(D,'minimize');
            D0.hyperparam = min_hyperparam{i,j,k,q};
            D = D0;
          end
          if ~isempty(min_hyperparam{i,j,k,q})
            D.hyperparam = min_hyperparam{i,j,k,q};
          end
        end
        
        % estimate lateralization index LI
        if calc_LI
          D.hemi = 'lh';
          if isfield(D,'eqdist')
            [BrainAGE_lh, ~, D] = cg_BrainAGE_GPR(D);
          else
            BrainAGE_lh = cg_BrainAGE_GPR(D);
          end

          D.hemi = 'rh';
          if isfield(D,'eqdist')
            [BrainAGE_rh, ~, D] = cg_BrainAGE_GPR(D);
          else
            BrainAGE_rh = cg_BrainAGE_GPR(D);
          end

          BrainAGE = (BrainAGE_lh - BrainAGE_rh)./(BrainAGE_lh + BrainAGE_rh);
        else
          if isfield(D,'eqdist')
            [BrainAGE, ~, D] = cg_BrainAGE_GPR(D);
          else
            BrainAGE = cg_BrainAGE_GPR(D);
          end
        end

        % move on if training failed
        if all(isnan(BrainAGE)) || std(BrainAGE)==0
          BrainAGE_all = BrainAGE;
          BA_unsorted = [BA_unsorted, BrainAGE];
          if isfield(D,'define_cov')
            PredictedAge_unsorted = [PredictedAge_unsorted, BrainAGE];
          else
            PredictedAge_unsorted = [PredictedAge_unsorted, BrainAGE + D.age_test];
          end
          % prepare groupwise data
          ind_groups = [];
          for o = 1:n_groups
            ind_groups = [ind_groups; D.ind_groups{o}];
          end
          
          BA = [BA, BrainAGE(ind_groups)];
          if isfield(D,'define_cov')
            EA = [EA, BrainAGE(ind_groups)];
          else
            EA = [EA, BrainAGE(ind_groups) + D.age_test(ind_groups)];
          end
          continue
        end
        
        % dont' apply trend correction during k-fold validation
        if D.trend_degree >= 0 && ~isfield(D,'k_fold')
          BrainAGE = apply_trend_correction(BrainAGE,D.age_test,D);
        end
        if isfield(D,'define_cov')
          PredictedAge = BrainAGE;
        else
          PredictedAge = BrainAGE + D.age_test;
        end

        RMSE = sqrt(sum((BrainAGE(D.ind_adjust)).^2)/length(D.ind_adjust));
        MAE  = mean(abs(BrainAGE(D.ind_adjust)));
        
        % correlation coefficient is not meaningful for too small samples
        if length(D.ind_adjust) > 5
          cc = corrcoef(PredictedAge(D.ind_adjust),D.age_test(D.ind_adjust));
        else
          cc(1,2) = NaN;
        end
        
        if isfield(D,'run_kfold') && D.run_kfold > 0
          if isfield(D,'run_repetition') && D.k_fold_reps > 1
            fprintf('repetition %02d/%02d fold %02d/%02d %s %s %s | test data: n=%d corr: %1.3f, MAE: %2.3f, RMSE: %2.3f\n',...
               D.run_repetition, D.k_fold_reps, D.run_kfold,D.k_fold,D.res,D.smooth,seg,numel(D.ind_adjust),cc(1,2),MAE,RMSE);
          else
            fprintf('fold %02d/%02d %s %s %s | test data: n=%d corr: %1.3f, MAE: %2.3f, RMSE: %2.3f\n',...
               D.run_kfold,D.k_fold,D.res,D.smooth,seg,numel(D.ind_adjust),cc(1,2),MAE,RMSE);
          end
        else
          fprintf('%s %s %s | test data: n=%d corr: %1.3f, MAE: %2.3f, RMSE: %2.3f\n',...
             D.res,D.smooth,seg,numel(D.ind_adjust),cc(1,2),MAE,RMSE);
        end
        
        data_cell = cell(1,n_groups);
        data = [];
        avg_BrainAGE = zeros(1,n_groups);
        age_test = zeros(n_groups,2);
        median_BrainAGE = zeros(1,n_groups);
        SD_BrainAGE = zeros(1,n_groups);
    
        % prepare group-wise data
        ind_groups = [];
        allow_violin = 1;
        for o = 1:n_groups
          if length(D.ind_groups{o}) < 10
            allow_violin = 0;
          end
          ind_groups = [ind_groups; D.ind_groups{o}];
          data_cell{o} = BrainAGE(D.ind_groups{o});
          avg_BrainAGE(o) = mean(BrainAGE(D.ind_groups{o}));
          age_test(o,1) = mean(D.age_test(D.ind_groups{o}));
          age_test(o,2) = std(D.age_test(D.ind_groups{o}));
          median_BrainAGE(o) = median(BrainAGE(D.ind_groups{o}));
          SD_BrainAGE(o) = std(BrainAGE(D.ind_groups{o}));
          data = [data; BrainAGE(D.ind_groups{o})];
        end
    
        % save data for all spatial resolutions and smoothing sizes
        BA_unsorted = [BA_unsorted, BrainAGE];
        if isfield(D,'define_cov')
          PredictedAge_unsorted = [PredictedAge_unsorted, BrainAGE];
          EA          = [EA, BrainAGE(ind_groups)];
        else
          PredictedAge_unsorted = [PredictedAge_unsorted, BrainAGE + D.age_test];
          EA          = [EA, BrainAGE(ind_groups) + D.age_test(ind_groups)];
        end
        BA          = [BA, data];
        
        if style == 1
          opt = struct('violin',2,'showdata',1,'groupcolor',groupcolor);
        else
          opt = struct('style',3,'groupcolor',groupcolor);
        end
        
        if ~allow_violin
          opt = struct('style',0,'groupcolor',groupcolor,'showdata',1);
          style = 1;
        end
        
        if D.verbose > 1
          figure(21)
          cat_plot_boxplot(data_cell,opt);
    
          set(gcf,'Name',name,'MenuBar','none');

          if style == 1
            set(gca,'XTick',1:n_groups,'XTickLabel',D.name_groups);
          else
            set(gca,'YTick',1:n_groups,'YTickLabel',D.name_groups(n_groups:-1:1,:));
          end
        end

        % print age of groups
        if D.verbose, fprintf('Age [years]:\n'); end
        if D.verbose
          fprintf('%20s\t','Group');
          for o = 1:n_groups
            fprintf('%20s\t',deblank(D.name_groups(o,:)));
          end
        
          fprintf('\n'); fprintf('%20s\t','Mean'); for o = 1:n_groups, fprintf('%20.3f\t',age_test(o,1)); end
          fprintf('\n'); fprintf('%20s\t','SD');   for o = 1:n_groups, fprintf('%20.3f\t',age_test(o,2)); end
          fprintf('\n%20s\t','Sample size');

          for o = 1:n_groups
            fprintf('%20s\t',sprintf('n=%d',length(D.ind_groups{o})));
          end
          fprintf('\n\n');

        end

        if D.verbose > 1
          if calc_LI
            % limit y-axis for LI, because there might be huge values sometimes
            if style == 1
              ylim([-5 5]);
              ylabel('Lateralization index of BrainAGE');
            else
              xlim([-5 5]);
              xlabel('Lateralization index of BrainAGE');
            end
            if D.verbose, fprintf('Lateralization index of BrainAGE:\n'); end
          else
            if style == 1
              ylabel('BrainAGE [years]');
            else
              xlabel('BrainAGE [years]');
            end
            if D.verbose, fprintf('BrainAGE [years]:\n'); end
          end
          set(gca,'FontSize',20);
        end
        
        % print BrainAGE of groups
        if D.verbose
          fprintf('%20s\t','Group');
          for o = 1:n_groups
            fprintf('%20s\t',deblank(D.name_groups(o,:)));
          end
        
          fprintf('\n'); fprintf('%20s\t','Mean');   for o = 1:n_groups, fprintf('%20.3f\t',avg_BrainAGE(o)); end
          fprintf('\n'); fprintf('%20s\t','Median'); for o = 1:n_groups, fprintf('%20.3f\t',median_BrainAGE(o)); end
          fprintf('\n'); fprintf('%20s\t','SD'); for o = 1:n_groups, fprintf('%20.3f\t',SD_BrainAGE(o)); end
          fprintf('\n');
    
        end 
  
        % ANOVA + T-Test
        if n_groups > 1
          if exist('cg_anova1')
            group = ones(length(D.ind_groups{1}),1);
            for o = 2:n_groups
              group = [group; o*ones(length(D.ind_groups{o}),1)];
            end
            Panova = cg_anova1(data,group,'off');
            if D.verbose, fprintf('\nANOVA P-value: %f\n',Panova); end
  
            P = zeros(n_groups,n_groups);
  
            for o = 1:n_groups
              for p = 1:n_groups
                [H,P(o,p)] = cg_ttest2(BrainAGE(D.ind_groups{o}),BrainAGE(D.ind_groups{p}),0.05,'left');
              end
            end
              
            fprintf('T-test P-value (one-tailed):\n');
            fprintf('%20s\t','Group');
            for o = 1:n_groups
              fprintf('%20s\t',deblank(D.name_groups(o,:)));
            end
            fprintf('\n');
  
            for o = 1:n_groups
              fprintf('%20s\t',deblank(D.name_groups(o,:)));
              for p = 1:n_groups
                if P(o,p) <= 0.001 || P(o,p) >= 0.999
                  fprintf('%20.7f***\t',P(o,p));
                elseif P(o,p) <= 0.01 || P(o,p) >= 0.99
                  fprintf('%20.7f **\t',P(o,p));
                elseif P(o,p) <= 0.05 || P(o,p) >= 0.95
                  fprintf('%20.7f  *\t',P(o,p));
                else
                  fprintf('%20.7f\t',P(o,p));
                end
              end
              fprintf('\n');
            end
        
            if ~isempty(find(P<=0.05))
              fprintf('****************************\n');
              fprintf('Significant result found\n');
              fprintf('****************************\n\n');
            end
                
          else
            fprintf('Warning: cg_anova1 not found.\n');
          end
  
        elseif isfield(D,'corr') % estimate correlation for one group
          [R, P] = corrcoef(BrainAGE(~isnan(D.corr)),D.corr(~isnan(D.corr)));
          fprintf('Correlation r=%g (P=%g)\n',R(1,2),P(1,2));
        end
        
      end
    end
  end        
end

% estimate ensembles
if multiple_BA && ((isfield(D,'run_kfold') && ~D.run_kfold) || ~isfield(D,'run_kfold'))

  if 0
    % apply initial trend correction to all models separately if not k-fold validation
    if D.trend_degree >= 0 && ~isfield(D,'k_fold')
      for i=1:size(BA_unsorted,2)
        BA_unsorted(:,i) = apply_trend_correction(BA_unsorted(:,i),D.age_test,D,0);
      end
    end
  end
  
  if D.ensemble_method == 4, fprintf('GPR ensemble approach is only useful for k-fold validation.\n'); end
  BA_unsorted_weighted  = ensemble_models(BA_unsorted,D.age_test,D);
    
  % apply final trend correction to weighted model if not k-fold validation
  if D.trend_degree >= 0 && ~isfield(D,'k_fold')
    BA_unsorted_weighted = apply_trend_correction(BA_unsorted_weighted,age,D);
  end
  
  BA_weighted = BA_unsorted_weighted(ind_groups);
  MAE_ctl_weighted = mean(abs(BA_unsorted_weighted(D.ind_adjust)));

  fprintf('-----------------------------------------------------------------\n'); 
  fprintf('weighted MAE (testdata): %g\n',MAE_ctl_weighted); 
  fprintf('-----------------------------------------------------------------\n'); 
  
  if D.verbose
    figure(22)
    plot(D.age_test(D.ind_adjust),D.age_test(D.ind_adjust)+BA_unsorted_weighted(D.ind_adjust),'.')
    hold on 
    line([0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],[0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],...
      'Color',[0 0 0],'lineWidth',2);
    hold off
    title('Test Data')
    xlabel('Chronological Age [years]')
    ylabel('Predicted Age [years]')

    figure(23)
    clf
    hold on 
    for i = 1:n_groups
      plot(D.age_test(D.ind_groups{i}),D.age_test(D.ind_groups{i})+BA_unsorted_weighted(D.ind_groups{i}),'*','color',groupcolor(i,:))
    end
    line([0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],[0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],...
      'Color',[0 0 0],'lineWidth',2);
    hold off
    title('Data')
    xlabel('Chronological Age [years]')
    ylabel('Predicted Age [years]')
    legend(D.name_groups,'Location','NorthWest')
    set(gca,'FontSize',20);

  end
  
  if n_groups > 1
    if exist('cg_anova1')
      Panova = cg_anova1(BA_weighted,group,'off');
      if D.verbose, fprintf('\nWeighted ANOVA P-value: %f\n',Panova); end
    
      P = zeros(n_groups,n_groups);
      for o = 1:n_groups
        for p = 1:n_groups
          [H,P(o,p)] = cg_ttest2(BA_unsorted_weighted(D.ind_groups{o}),BA_unsorted_weighted(D.ind_groups{p}),0.05,'left');
        end
      end
    
      fprintf('T-test P-value (one-tailed):\n');
      fprintf('%20s\t','Group');
      for o = 1:n_groups
        fprintf('%20s\t',deblank(D.name_groups(o,:)));
      end
      fprintf('\n');
    
      for o = 1:n_groups
        fprintf('%20s\t',deblank(D.name_groups(o,:)));
        for p = 1:n_groups
          if P(o,p) <= 0.001 || P(o,p) >= 0.999
            fprintf('%20.7f***\t',P(o,p));
          elseif P(o,p) <= 0.01 || P(o,p) >= 0.99
            fprintf('%20.7f **\t',P(o,p));
          elseif P(o,p) <= 0.05 || P(o,p) >= 0.95
            fprintf('%20.7f  *\t',P(o,p));
          else
            fprintf('%20.7f\t',P(o,p));
          end
        end
        fprintf('\n');
      end
    
      if ~isempty(find(P<=0.05))
        fprintf('****************************\n');
        fprintf('Significant result found\n');
        fprintf('****************************\n\n');
      end
  
    else
      fprintf('cg_anova1 not found.\n');
    end
  end
  BrainAGE = BA_weighted;
  BrainAGE_unsorted = BA_unsorted_weighted;
 
  if nargout > 2
    if isfield(D,'define_cov')
      BrainAGE_all = PredictedAge_unsorted;
    else
      BrainAGE_all = PredictedAge_unsorted - D.age_test;
    end
    BrainAGE_all = BrainAGE_all(ind_groups,:);
  end
else
  BrainAGE = BA;
  BrainAGE_unsorted = BA_unsorted;
  if nargout > 2
    BrainAGE_all = BA;
  end
end

% show plot for multiple values if defined
if multiple_BA && ((isfield(D,'run_kfold') && ~D.run_kfold) || ~isfield(D,'run_kfold'))
      
  ind_groups = [];
  for o = 1:n_groups
    ind_groups = [ind_groups; D.ind_groups{o}];
    data_cell{o} = BA_unsorted_weighted(D.ind_groups{o});
    avg_BrainAGE(o) = mean(BA_unsorted_weighted(D.ind_groups{o}));
    age_test(o,1) = mean(D.age_test(D.ind_groups{o}));
    age_test(o,2) = std(D.age_test(D.ind_groups{o}));
    median_BrainAGE(o) = median(BA_unsorted_weighted(D.ind_groups{o}));
  end
    
  % print BrainAGE of groups
  if D.verbose && all(sum(isnan(BA_unsorted)) == 0)
    fprintf('%20s\t','Group');
    for o = 1:n_groups
      fprintf('%20s\t',deblank(D.name_groups(o,:)));
    end
  
    fprintf('\n'); fprintf('%20s\t','Mean');   for o = 1:n_groups, fprintf('%20.3f\t',avg_BrainAGE(o)); end
    fprintf('\n'); fprintf('%20s\t','Median'); for o = 1:n_groups, fprintf('%20.3f\t',median_BrainAGE(o)); end
    fprintf('\n');

    figure(24)
    cat_plot_boxplot(data_cell,opt);

    set(gcf,'Name','Weighted BrainAGE','MenuBar','none');

    if style == 1
      set(gca,'XTick',1:n_groups,'XTickLabel',D.name_groups);
    else
      set(gca,'YTick',1:n_groups,'YTickLabel',D.name_groups(n_groups:-1:1,:));
    end
    set(gca,'FontSize',20);
    if style == 1
      ylabel('BrainAGE [years]');
    else
      xlabel('BrainAGE [years]');
    end
  end
end


%-------------------------------------------------------------------------------
function [BA_weighted, PredictedAge_weighted] = ensemble_models(BA,age_all,D,ind_test_array,ind_train_array)
%-------------------------------------------------------------------------------
% Estimate ensembles to combine different models
% and calculate weighted BrainAGE score
% Weighting is only estimated using data of control group that is indicated by 
% D.ind_adjust (i.e. the first group if not otherwise defined)

if size(BA,1) ~= size(age_all,1)
  age = age_all(D.ind_adjust);
else
  age = age_all;
end

% only estimates weights if BA is given from different models
if size(BA,2) == 1
  BA_weighted = BA;
  if isfield(D,'define_cov')
    PredictedAge_weighted = BA_weighted;
  else
    PredictedAge_weighted = BA_weighted + age;
  end
  fprintf('Disable model ensemble because we just have one model.\n');
  return
end

% use median as default if no ensemble method is given
if isfield(D,'ensemble')
  ensemble_method = D.ensemble;
else
  ensemble_method = 2;
end

if nargin < 4 && ensemble_method == 4
  fprintf('For GPR ensemble method you have to define 5 arguments.\n');
  ensemble_method = 1;
end

% remove columns with NaNs where model estimation failed
BA_corrected = BA;
[indx,indy] = find(isnan(BA_corrected));
if ~isempty(indy)
  indy = unique(indy);
  BA_corrected(:,indy) = [];
end

if isfield(D,'define_cov')
  PredictedAge_corrected = BA_corrected;
else
  PredictedAge_corrected = BA_corrected + age;
end

switch ensemble_method
case 0   % use model with lowest MAE
  
  BA_ind = BA_corrected(D.ind_adjust,:);
  [mn, mn_ind] = min(mean(abs(BA_ind)));
  
  fprintf('\nModel with lowest MAE is %d\n',mn_ind);
  
  weighted_PredictedAge = PredictedAge_corrected(:,mn_ind);
  BA_weighted = weighted_PredictedAge - age;

case 1   % use GLM estimation to minimize MAE
  
  PredictedAge_ind = PredictedAge_corrected(D.ind_adjust,:);
  Beta = pinv(PredictedAge_ind)*age(D.ind_adjust);
  
  fprintf('\nEstimated Betas using %d subjects:\t',numel(D.ind_adjust));
  fprintf('%.2f ',Beta);
  fprintf('\n');
  
  weighted_PredictedAge = PredictedAge_corrected*Beta;
  BA_weighted = weighted_PredictedAge - age;

case 2 % simple median of all models
  
  BA_weighted = median(BA_corrected,2);
case 3   % use GLM estimation to maximize group differences or correlation (EXPERIMENTAL!)

  if ~isfield(D,'contrast')
    error('D.contrast has to be defined.');
  end

  % F-contrast ?
  if size(D.contrast,1) == numel(D.ind_groups) && size(D.contrast,2) == numel(D.ind_groups)
    c = D.contrast;
    D.contrast = [];
    ind = [];
    for i=1:numel(D.ind_groups)
      D.contrast = [D.contrast; repmat(c(i,:),numel(D.ind_groups{i}),1)];
      ind = [ind; D.ind_groups{i}];
    end
    % re-order w.r.t. D.ind_groups
    if max(ind) == numel(ind)
      D.contrast = D.contrast(ind,:);
    else
      fprintf('Warning: Order of D.ind_groups should be not mixed because we cannot identify the correct order!\n');
    end
    group_diff = 0;
  elseif size(D.contrast,1) == size(BA_corrected,1)
    group_diff = 0;
  elseif size(D.contrast,2) == numel(D.ind_groups)
    group_diff = 1;
  else
    error('D.contrast has different size than number of groups.');    
  end
  
  warning('Experimental: Not yet working properly!');

  if group_diff
    X = []; Y = [];
    for i=1:size(D.contrast,2)
      if any(D.contrast(:,i) ~= 0)
        X = [X; D.contrast(:,i)*ones(numel(D.ind_groups{i}),1)];  
        Y = [Y; BA_corrected(D.ind_groups{i},:)];  
      end
    end
  else
    X = D.contrast;
    Y = BA_corrected;
  end
  
  % use only those data that are defined in ind_groups
  ind = [];
  for i=1:numel(D.ind_groups)
    ind = [ind; D.ind_groups{i}];
  end
  
  if group_diff
    Beta = pinv(Y)*X;
  else
    % we have to excluded NaNs for Beta estimation
    ind_finite = isfinite(X);
    Beta = pinv(Y(ind(ind_finite),:))*X(ind_finite);
  end
  Beta = Beta./sum(Beta); % ensure that sum is 1

  fprintf('\nPredicted Betas using %d subjects:\t',numel(D.ind_adjust));
  fprintf('%.2f ',Beta);
  fprintf('\n');

  BA_weighted = sum(BA_corrected*Beta,2);

  % scale estimated BA values by ratio between SD of original and estimated BA values to get the same range
  BA_weighted = BA_weighted*mean(std(BA_corrected))/mean(std(BA_weighted));
case 4   % Stacking: use GPR to combine models
        
  PredictedAge_ind = PredictedAge_corrected;
  
  hyperparam = struct('mean', 5, 'cov', [5 4 -2], 'lik', 1);
%  hyperparam = minimize(hyperparam, @gp, -100, @infGaussLik, @meanConst, @covRQiso, @likGauss, PredictedAge_ind, age)
  for i=1:numel(ind_test_array)
    
    % get indices for training and test and only include data that are
    % defined in D.ind_adjust
    ind_train = ind_train_array{i};
    ind_train = ind_train(ismember(ind_train,D.ind_adjust));
    ind_test  = ind_test_array{i};
    ind_test  = ind_test(ismember(ind_test,D.ind_adjust));
    
    Y_train = PredictedAge_ind(ind_train,:);
    Y_test  = PredictedAge_ind(ind_test,:);
    
    % Regression using GPR
    % yfit is the predicted mean, and ys is the predicted variance
    PredictedAge = gp(hyperparam, @infGaussLik, @meanConst, @covRQiso, @likGauss, Y_train, age(ind_train), Y_test);
    BA_weighted(ind_test) = PredictedAge-age(ind_test);
    
  end
  
  BA_weighted = BA_weighted';
case 5 % weighted average of all models w.r.t. MAE
  MAE = mean(abs(BA_corrected));
  weight = 1./MAE.^2;
  weight = weight/mean(weight);
  PredictedAge_weighted = mean(PredictedAge_corrected.*weight,2);
  BA_weighted = PredictedAge_weighted - age;
end

if isfield(D,'define_cov')
  PredictedAge_weighted = BA_weighted;
else
  PredictedAge_weighted = BA_weighted + age;
end

%-------------------------------------------------------------------------------
function [BrainAGE, PredictedAge, Adjustment] = apply_trend_correction(BrainAGE,age,D,verbose)
%-------------------------------------------------------------------------------
% BrainAGE     - adjusted BrainAGE
% PredictedAge - adjusted PredictedAge
% Adjustment   - array of adjustments for indexed data

PredictedAge = BrainAGE + age;
if D.trend_degree < 0
  return
end

% be verbose by default
if nargin < 4
  verbose = 1;
end

% if site_adjust is empty we don't apply site adjustment
if isempty(D.site_adjust)
  site_adjust = ones(D.n_data,1);
else
  site_adjust = D.site_adjust;
end

Adjustment = cell(max(site_adjust),1);

% save MAE before correction
MAE0 = mean(abs(BrainAGE));

% apply trend correction for each site separately
for i = 1:max(site_adjust)
  ind_site = find(site_adjust == i);
  
  % build index for each site
  ind_site_adjust = [];
  for j=1:length(D.ind_adjust)
    if site_adjust(D.ind_adjust(j)) == i
      ind_site_adjust = [ind_site_adjust; D.ind_adjust(j)];
    end
  end

  % apply trend correction for indexed data
  if ~isempty(ind_site_adjust)
  
    % test whether sample size is too small
    if length(ind_site_adjust) < 15
      if verbose, fprintf('Sample #%d for site-specific trend-correction is too small (n=%d). Use only offset-correction by Median instead.\n',...
          i,length(ind_site_adjust)); end
      offset = median(BrainAGE(ind_site_adjust));
      BrainAGE(ind_site) = BrainAGE(ind_site) - offset;
      PredictedAge(ind_site) = PredictedAge(ind_site) - offset;
      
    % test whether age range is too small
    elseif (max(age(ind_site_adjust)) - min(age(ind_site_adjust))) < 1
      if verbose, fprintf('Age range for site-specific trend-correction is too small (%g). Use only offset-correction by Median instead.\n',...
          (max(age(ind_site_adjust)) - min(age(ind_site_adjust)))); end
      offset = median(BrainAGE(ind_site_adjust));
      BrainAGE(ind_site) = BrainAGE(ind_site) - offset;
      PredictedAge(ind_site) = PredictedAge(ind_site) - offset;
      
    % everything else
    else

      % use BrainAGE for obtaining trend
      if D.trend_method == 1
        Y = BrainAGE;
      else
        % use predicted age for obtaining trend
        Y = PredictedAge;
      end
      
      % polynomial trend
      if D.trend_degree > 0
        G = [ones(length(age(ind_site)),1) cg_polynomial(age(ind_site),D.trend_degree)];
        G_indexed = [ones(length(age(ind_site_adjust)),1) cg_polynomial(age(ind_site_adjust),D.trend_degree)];
      % use just offset
      else
        G = ones(length(age(ind_site)),1);
        G_indexed = ones(length(age(ind_site_adjust)),1);
      end
      if verbose, fprintf('Remove trend degree %d using %d subjects of site %d.\n',D.trend_degree,length(ind_site_adjust),i); end
        
      % estimate beta only for indexed data (e.g. control subjects)
      Beta = pinv(G_indexed)*Y(ind_site_adjust);
      
      % and remove effects for all data
      GBeta = G*Beta;

      BrainAGE0 = BrainAGE;
      
      % use BrainAGE for obtaining trend
      if D.trend_method == 1
        PredictedAge(ind_site)  = PredictedAge(ind_site)  - GBeta;
        BrainAGE = PredictedAge - age;
        Adjustment{i} = GBeta;
      else
        % use predicted age for obtaining trend
        PredictedAge(ind_site)  = (PredictedAge(ind_site) - Beta(1));
        if D.trend_degree == 1
          PredictedAge(ind_site) = PredictedAge(ind_site)/Beta(2);
        end
        BrainAGE = PredictedAge - age;
      end
      Adjustment{i} = BrainAGE0(ind_site) - BrainAGE(ind_site);
    end
  else
    if D.trend_degree >= 0
      if verbose, fprintf('Warning: No subjects found in sample #%d for site-specific trend-correction\n',i); end
    end
  end
end

avg_BrainAGE = mean(BrainAGE(D.ind_adjust));
BrainAGE = BrainAGE - avg_BrainAGE;
PredictedAge = PredictedAge - avg_BrainAGE;

MAE = mean(abs(BrainAGE));

if MAE0/MAE > 4
  if verbose, fprintf('Warning: Large discrepancy between MAE before and after correction which points to a too narrow age range of training data!\n'); end
end

function C = nejm
  C = [
    '#BC3C29'
    '#0072B5'
    '#E18727'
    '#20854E'
    '#7876B1'
    '#6F99AD'
    '#FFDC91'
    '#EE4C97'
    '#8C564B'
    '#BCBD22'
    '#00A1D5'
    '#374E55'
    '#003C67'
    '#8F7700'
    '#7F7F7F'    
    '#353535'    
  ];
  C = reshape(sscanf(C(:,2:end)','%2x'),3,[]).'/255;

