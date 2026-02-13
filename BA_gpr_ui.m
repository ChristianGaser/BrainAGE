function [BrainAGE, BrainAGE_unsorted, BrainAGE_all, D, age] = BA_gpr_ui(D)
% [BrainAGE, BrainAGE_unsorted, BrainAGE_all, D] = BA_gpr_ui(D)
% User interface for BrainAGE estimation in BA_gpr.m
%
% D.data            - test sample for BrainAGE estimation
% D.seg_array       - segmentations
%                     {'rp1'} use GM
%                     {'rp2'} use WM
%                     {'rp1,'rp2'} use GM or WM
%                     {'rp1+rp2'} use both GM+WM (by concatinating both)
% D.res_array       - spatial resolution of data as char or cell (can be a cell array to try different values: {'4','8'})
% D.smooth_array    - smoothing size as char or cell (can be a cell array to try different values: {'s4','s8'})
% D.train_array     - cell array name(s) of training samples
%     Healthy adults:
%            IXI547 - IXI
%          OASIS316 - OASIS
%       OASIS3_1752 - OASIS3 (all time points of 549 subjects)
%        OASIS3_549 - OASIS3 (only last time point)
%         CamCan651 - CamCan
%           SALD494 - SALD
%           NKIe629 - NKIe (minimum age 6y)
%           NKIe516 - NKIe (minimum age 18y)
%     ADNI231Normal - ADNI Normal-sc-1.5T 
%
%     Children data:
%            NIH394 - NIH objective 1
%            NIH879 - NIH release 4.0
%            NIH755 - NIH release 5.0
%
%     Children + adults:
%         fCONN772  - fcon-1000 (8-85 years)
%
% D.relnumber       - CAT12 release (e.g. '_CAT12.9')
% D.age_range       - age range of training data
%                     This is a very important and useful parameter that should be adapted to your data.
%                     If your sample consists mainly of elderly subjects or patients with dementia, 
%                     it is recommended to use a higher minimum age range for the training data to 
%                     account for aging effects in subjects older than 65 years, where the largest 
%                     aging effects occur (e.g. [50 Inf]). 
%                     In the case of much younger subjects (20..35 years), an age range of 18..50 
%                     (or even 60) may be more appropriate to account for the different age trajectories 
%                     at younger ages.
%                     [0 Inf] use all data (default)
%                     [50 80] use age range of 50..80
% D.ind_adjust      - define indices for adjusting data according to trend defined with D.trend_degree
%                     usually this is the control group that is used in order to adjust the data
% D.site_adjust     - If data are acquired at different sites (e.g. using different scanners or sequences) the trend 
%                     correction should be estimated and applied for each site seperately. A vector with coding of the scanners is required.
%                     If this parameter is empty this ensures that no site-specific adjustment is made even for multiple
%                     training data that are combined with a "+". 
% D.ind_groups      - define indices of groups, if D.ind_adjust is not given, then the first group of this index will
%                     be used for adjusting data according to trend defined with D.trend_degree
% D.trend_method    - use different methods for estimating trend:
%                     0 skip trend correction (set trend_degree to -1)
%                     1 use BrainAGE for trend correction (as used in Smith et al. 2019, default)
%                     2 use predicted age for trend correction (as used in Cole et al. 2018)
% D.trend_degree    - estimate trend with defined order using healthy controls and apply it to all data (set to -1 for skipping trend correction)
% D.trend_ensemble  - apply trend correction for each ensemble separately before bagging/stacking
%                     0 skip trend correction for each ensemble (default)
%                     1 apply trend correction for each ensemble
% D.hyperparam      - GPR hyperparameters (.mean and .lik)
% D.RVR             - use old RVR method
% D.PCA             - apply PCA as feature reduction (default=1), values > 1 define number of PCA components
% D.PCA_method      - method for PCA
%                    'eig' Eigenvalue Decomposition of the covariance matrix (faster but less accurate, for compatibiliy)
%                    'svd' Singular Value Decomposition of X (the default)
% D.p_dropout       - Dropout probability to randomly exclude voxels/data points to implement an uncertainty-aware approach using a 
%                     Monte-Carlo Dropout during inference. That means that during testing, voxels are randomly dropped out according 
%                     to the dropout probabilities. This process is repeated multiple times, and each time, the model produces 
%                     a different output. By averaging these outputs, we can obtain a more robust prediction and estimate the model's 
%                     uncertainty in its predictions. A meaningful dropout probability is 0.1, which means that 10% of the data points 
%                     are excluded. The default is 0.
% D.ensemble        - ensemble method to combine different models
%                     0 - Majority voting: use model with lowest MAE
%                     1 - Weighted GLM average: use GLM estimation to estimate model weights to minimize MAE
%                     2 - Average: use mean to weight different models
%                     3 - GLM: use GLM estimation to maximize variance to a group or a regression parameter (EXPERIMENTAL!)
%                     5 - Weighted Average: (average models with weighting w.r.t. squared MAE) (default)
%                     6 - GLM: use GLM estimation for different tissues (i.e. GM/WM) to maximize variance to a group or a regression parameter (EXPERIMENTAL!)
%                         In contrast to ensemble model 3, we here only use the mean tissue values and not all models to estimate weights
% D.contrast        - define contrast to maximize group differences (use only if D.ensemble is 3 or 6) (e.g. [1 -1])
%                     D.contrast can be also a vector which is used to maximize variance between BrainAGE and this parameter.
% D.dir             - directory for databases and code
% D.verbose         - verbose level
%                     0 - suppress long outputs
%                     1 - print meaningful outputs (default)
%                     2 - print long outputs
% D.parcellation    - use parcellation into lobes to additionally estimate local BrainAGE values:
%                     https://figshare.com/articles/dataset/Brain_Lobes_Atlas/971058
%                     0 - estimate global BrainAGE (default)
%                     1 - estimate local BrainAGE for different lobes for both hemispheres
% D.threshold_std   - all data with a standard deviation > D.threshold_std of mean covariance are excluded (after covarying out effects of age)
%                     meaningful values are 1,2 or Inf
% D.eqdist          - options for age and sex equalization between test and train
% D.eqdist.weight   - vector of size 2 that allows to weight the cost function for age and sex equalization
% D.eqdist.range    - matrix 2 x 2 which defines the age range and sex range for equalization
% D.eqdist.tol      - vector of size 2 that defines tolerance between mean value of age_test and age_train and male_test and male_train
%                     A recommended value is [3 Inf], which means that a mean difference in age of 3 years is allowd,
%                     while sex is not considered.
% D.eqdist.debug    - print debug info if set (default 0)
% D.corr            - additionally define parameter that can be correlated to BrainAGE if only one group is given
% D.define_cov      - optionally define continous parameter that should be used instead of age for more general use 
%                     of GPR not only limited to BrainAGE
% D.style           - plot-style: 1: old style with vertical violin-plot; 2: new style with horizontal density plot (default)
% D.groupcolor      - matrix with (group)-bar-color(s), use jet(numel(data)) or other color functions (nejm by default)
% D.normalize_BA    - normalize BA values w.r.t. MAE to make BA less dependent from training sample (i.e. size) and scale 
%                     it to MAE of 5
% D.comcat          - If data are acquired at different sites (e.g. using different scanners or sequences) we can
%                     harmonize data using ComCAT. A vector with coding of the scanners is required for the test data (EXPERIMENTAL!).
% D.nuisance        - additionally define nuisance parameter for covarying out (e.g. gender)
% D.spiderplot.func - show spider (radar) plot either with mean or median values (only valid if D.parcellation is used):
%                     'median' - use median values 
%                     'mean'   - use mean values (default)
% D.spiderplot.range- range for spiderplot (default automatically find range)
%
% Parameter search
% ---------------
% Some selected parameters can be also defined as ranges to try different parameter settings.
% Examples:
% D.trend_degree = 1;
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
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

global min_hyperparam

% add cat12 path if not already done
if ~exist('cat_stat_polynomial','file')
  addpath(fullfile(spm('dir'),'toolbox','cat12'));
end

% use normalized BA as default and scale it to MAE of 5
if ~isfield(D,'normalize_BA')
  D.normalize_BA = 0;
else
  if D.normalize_BA == 1
    D.normalize_BA = 5;
  end
end

% check whether k-fold validation is defined
if isfield(D,'n_fold') || isfield(D,'k_fold')
  fprintf('-------------------------------------------------------------\n');
  fprintf('For using k-fold validation you have to call BA_gpr_kfold_ui.\n');
  fprintf('-------------------------------------------------------------\n');
  [BrainAGE, BrainAGE_unsorted, BrainAGE_all, D, age] = BA_gpr_kfold_ui(D);
  return;
end

if ~isfield(D,'spiderplot') || (isfield(D,'spiderplot') && ~isfield(D.spiderplot,'func'))
  D.spiderplot.func = 'mean';
end

if ~isfield(D,'trend_degree')
  D.trend_degree = 1;
end

if ~isfield(D,'trend_method')
  D.trend_method = 1;
end

if ~isfield(D,'trend_ensemble')
  D.trend_ensemble = 0;
end

if ~isfield(D,'ensemble')
  D.ensemble = 5;
end

if D.ensemble < 0 && ~exist('fmincon')
  fprintf('In order to use non-linear optimization you need the Optimization Toolbox.\n');
  return
end

if ~isfield(D,'age_range')
  D.age_range = [0 Inf];
end

if ~isfield(D,'hyperparam')
  D.hyperparam = struct('mean', 100, 'lik', -1);
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

if ~isfield(D,'RVR')
  D.RVR = 0;
end

if ~isfield(D,'PCA_method')
  D.PCA_method = 'svd';
end

if ~isfield(D,'parcellation')
  D.parcellation = 0;
end

if D.trend_method > 1 && D.trend_degree > 1
  D.trend_degree = 1;
  fprintf('Only use linear trend correction for method that uses predicted age for obtaining trend correction.\n');
end

if ~D.trend_method
  if D.trend_degree > -1, fprintf('Disable trend correction.\n'); end
  D.trend_degree = -1;
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

% set default for droput probability 
if ~isfield(D,'p_dropout')
  D.p_dropout = 0;
end

if iscell(D.data)
  D.data = char(D.data);
end

% if comcat is defined and set to 0 then remove this field
if isfield(D,'comcat') && isscalar(D.comcat) &&  D.comcat == 0
  D = rmfield(D,'comcat');
end

region_names = {'R Frontal','R Parietal','R Occipital','R Temporal','R Subcortical/Cerebellum',...
         'L Frontal','L Parietal','L Occipital','L Temporal','L Subcortical/Cerebellum'};

region_names_surf = {'R Frontal','R Parietal','R Occipital','R Temporal',...
         'L Frontal','L Parietal','L Occipital','L Temporal'};

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
    % data size from that
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
    
    for l = 1:n_train
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

if isfield(D,'weighting') && ~isfield(D,'ensemble')
  D.ensemble = D.weighting;
  fprintf('This option is deprecated. Use the option ''ensemble'' instead.\n');
  return
end

% check whether contrast was defined for ensemble=3
if isfield(D,'ensemble')
  if (abs(D.ensemble(1)) == 3 || abs(D.ensemble(1)) == 6) && ~isfield(D,'contrast')
    error('D.contrast has to be defined.');
  end
end

% print some parameters
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
if D.p_dropout
  fprintf('Prob-Dropout: \t%d\n',D.p_dropout);
end
if D.RVR
  fprintf('RVR:          \t%d\n',D.RVR);
end
fprintf('PCA:           \t%d (method: %s)\n',D.PCA,D.PCA_method);
fprintf('Trend method:  \t%d\n',D.trend_method);
fprintf('Age-Range:     \t%g-%g\n',D.age_range(1),D.age_range(2));
fprintf('--------------------------------------------------------------\n');
if isfield(D,'parcellation') && D.parcellation
  fprintf('Estimate local BrainAGE with parcellation into lobes.\n');
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
          for l = 1:n_train
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
          for l = 1:n_seg
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
          for l = 1:n_train
            name = [D.smooth D.seg{1} '_' D.res 'mm_' D.data(ind_plus(l)+1:ind_plus(l+1)-1) D.relnumber];
            if D.verbose > 1, fprintf('BA_gpr_ui: load %s\n',name); end
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
          if contains(D.seg{1}, 'mesh')
            name = [D.smooth '.' D.seg{1} '_' D.data D.relnumber];
          else
            name = [D.smooth D.seg{1} '_' D.res 'mm_' D.data D.relnumber];
          end
          if strcmp(name(1),'.')
            name = name(2:end);
          end
          if D.verbose > 1, fprintf('BA_gpr_ui: load %s\n',name); end
          load(name);
        end
        
        n_data = size(Y,1);
        if D.n_data ~= n_data
          fprintf('\n-----------------------------------------------------------------\n'); 
          fprintf('Data size differs for %s (%d vs. %d)',name,D.n_data,n_data);
          fprintf('\n-----------------------------------------------------------------\n'); 
          return
        end
    
        if D.verbose > 1, fprintf('\n-----------------------------------------------------------------\n%s\n',name); end
        D.Y_test = single(Y); clear Y
        if isfield(D,'define_cov')
          age = D.define_cov;
          male = ones(size(age));
        end
                
        D.age_test = age;

        if exist('male','var')
          D.male_test = male;
        end
    
        % use additional segmentation if defined
        for l = 2:n_seg
          % find potential "+" indicating to combine data
          ind_plus = strfind(D.data,'+');
          if ~isempty(ind_plus)
            Y0    = [];
            ind_plus = [0 ind_plus length(D.data)+1];
            n_train = numel(ind_plus)-1;
            for m = 1:n_train
              if contains(D.seg{l}, 'mesh')
                name = [D.smooth D.seg{l} '_' D.data(ind_plus(m)+1:ind_plus(m+1)-1) D.relnumber];
              else
                name = [D.smooth D.seg{l} '_' D.res 'mm_' D.data(ind_plus(m)+1:ind_plus(m+1)-1) D.relnumber];
              end
              if strcmp(name(1),'.')
                name = name(2:end);
              end
              if D.verbose > 1, fprintf('BA_gpr_ui: load %s\n',name); end
              load(name);
              Y0    = [Y0; single(Y)]; clear Y
            end
            D.Y_test = [D.Y_test Y0]; clear Y0
          else
            if contains(D.seg{l}, 'mesh')
              name = [D.smooth D.seg{l} '_' D.data D.relnumber];
            else
              name = [D.smooth D.seg{l} '_' D.res 'mm_' D.data D.relnumber];
            end
            if strcmp(name(1),'.')
              name = name(2:end);
            end
            if D.verbose > 1, fprintf('BA_gpr_ui: load %s\n',name); end
            load(name);
            D.Y_test = [D.Y_test single(Y)]; clear Y
          end
        end

        if D.verbose > 1, fprintf('\n'); end
                
        if ~isfield(D,'ind_groups')
          D.ind_groups = {1:length(D.age_test)};
        end
        
        n_groups = numel(D.ind_groups);
        
        if ~isfield(D,'groupcolor')
          groupcolor = cat_io_colormaps('nejm',n_groups);
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
                                

        if numel(D.age_range) ~=2
          error('Age range has to be defined by two values (min/max)');
        end
        
        % add spatial resolution to atlas name
        if isfield(D,'parcellation') && D.parcellation
          is_surf = contains(seg,'mesh');
          if is_surf
            lh_atlas = 'lh.Brain_Lobes.annot';
            [~, lh_label, ~] = cat_io_FreeSurfer('read_annotation',lh_atlas);
            rh_atlas = 'rh.Brain_Lobes.annot';
            [~, rh_label, ~] = cat_io_FreeSurfer('read_annotation',rh_atlas);

            % since the atlas regions are equal and go up to 200 we have to
            % multiply rh with 5
            atlas = [5*rh_label; lh_label];

            % lobes have higher numbers and start with 50
            regions = unique(atlas(atlas >= 50));

            % subcortical regions and cerebellum have to be removed (200
            % and 200*5)
            regions(regions==200)  = [];
            regions(regions==1000) = [];
          else
            atlas_name = ['Brain_Lobes_' D.res 'mm.mat'];          
            load(atlas_name)
            
            if ~exist('atlas','var')
              error('Atlas must contain atlas as variable');
            end
            regions = unique(atlas(atlas > 0));
          end
          
          D.n_regions = numel(regions);
          BrainAGE = [];
          for r = 1:D.n_regions
            if is_surf
              % for surface we have to use areas that are defined by ind
              D.mask = atlas(ind) == regions(r);
            else
              D.mask = atlas == regions(r);
            end
            [tmp, ~, D] = BA_gpr(D);
            BrainAGE = [BrainAGE tmp];
          end
        else
          [BrainAGE, ~, D] = BA_gpr(D);
          D.n_regions = 1;
        end
           
        % print information about training sample only once for D.threshold_std == Inf or otherwise always
        if D.verbose && (i==1 && j==1 && k==1) || isfinite(D.threshold_std)
          % only use selected data
          age_test = age(ind_test);
          fprintf('\n%d subjects used for training (age %3.1f..%3.1f years)\n',length(D.age_train),min(D.age_train),max(D.age_train));
          fprintf('Mean age\t%g (SD %g) years\nMales/Females\t%d/%d\n',mean(D.age_train),std(D.age_train),sum(D.male_train),length(D.age_train)-sum(D.male_train));
          if ~isfinite(D.threshold_std)
            fprintf('\n');
          end
          fprintf('\n%d subjects used for prediction (age %3.1f..%3.1f years)\n',length(age_test),min(age_test),max(age_test));
          if exist('male','var')
            fprintf('Mean age\t%g (SD %g) years\nMales/Females\t%d/%d\n',mean(age_test),std(age_test),sum(male),length(age_test)-sum(male));
          else
            fprintf('Mean age\t%g (SD %g) years\n',mean(age_test),std(age_test));
          end
          if ~isfinite(D.threshold_std)
            fprintf('\n');
          end

          figure(25)
          X0 = linspace(min([D.age_train; age_test]), max([D.age_train; age_test]), 20);
          H0 = histogram(D.age_train, X0);
          hold on
          H1 = histogram(age_test, X0);
          legend('Age training','Age test')
          title('Age distribution')
          set(gcf,'MenuBar','none');
          hold off
        end

        % move on if training failed for global BrainAGE
        if size(BrainAGE,2) == 1 && (all(isnan(BrainAGE)) || std(BrainAGE)==0)
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
        
        % dont' apply trend correction
        if (D.trend_ensemble || ~multiple_BA) && D.trend_degree >= 0
          BrainAGE = apply_trend_correction(BrainAGE,D.age_test,D);
        end
        if isfield(D,'define_cov')
          PredictedAge = BrainAGE;
        else
          PredictedAge = BrainAGE + D.age_test;
        end

        MAE  = mean(abs(BrainAGE(D.ind_adjust,:)));
        
        % correlation coefficient is not meaningful for too small samples
        % and regional estimates
        if length(D.ind_adjust) > 5
          cc = corrcoef([PredictedAge(D.ind_adjust,:) D.age_test(D.ind_adjust)]);
        else
          cc(1,2) = NaN;
        end
        
        fprintf('%s %s %s | test data: n=%d MAE/corr:',...
           D.res,D.smooth,seg,numel(D.ind_adjust));
        for r = 1:D.n_regions
          fprintf(' %2.3f/%1.3f',MAE(r),cc(end,r));
        end
        fprintf('\n');
        
        data_cell = cell(1,n_groups);
        data = [];
        avg_BrainAGE = zeros(n_groups,D.n_regions);
        age_test = zeros(n_groups,2);
        median_BrainAGE = zeros(n_groups,D.n_regions);
        SD_BrainAGE = zeros(n_groups,D.n_regions);
    
        % prepare group-wise data
        ind_groups = [];
        allow_violin = 1;
        for o = 1:n_groups
          if length(D.ind_groups{o}) < 10
            allow_violin = 0;
          end
          ind_groups = [ind_groups; D.ind_groups{o}];
          data_cell{o} = BrainAGE(D.ind_groups{o},:);
          avg_BrainAGE(o,:) = mean(BrainAGE(D.ind_groups{o},:));
          age_test(o,1) = mean(D.age_test(D.ind_groups{o},:));
          age_test(o,2) = std(D.age_test(D.ind_groups{o},:));
          median_BrainAGE(o,:) = median(BrainAGE(D.ind_groups{o},:));
          SD_BrainAGE(o,:) = std(BrainAGE(D.ind_groups{o},:));
          data = [data; BrainAGE(D.ind_groups{o},:)];
        end
    
        % save data for all spatial resolutions and smoothing sizes
        BA_unsorted = [BA_unsorted, BrainAGE];
        if isfield(D,'define_cov')
          PredictedAge_unsorted = [PredictedAge_unsorted, BrainAGE];
          EA = [EA, BrainAGE(ind_groups)];
        else
          PredictedAge_unsorted = [PredictedAge_unsorted, BrainAGE + D.age_test];
          EA = [EA, BrainAGE(ind_groups) + D.age_test(ind_groups)];
        end
        BA = [BA, data];
        
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
          set(gcf,'Name',name,'MenuBar','none');

          cat_plot_boxplot(data_cell,opt);
    
          if style == 1
            set(gca,'XTick',1:n_groups,'XTickLabel',D.name_groups);
          else
            set(gca,'YTick',1:n_groups,'YTickLabel',D.name_groups(n_groups:-1:1,:));
          end
        end

        % print age of groups
        if D.verbose > 1
          fprintf('Age [years]:\n');
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
          if style == 1
            ylabel('BrainAGE [years]');
          else
            xlabel('BrainAGE [years]');
          end
          set(gca,'FontSize',20);
        end
        
        % print BrainAGE of groups
        if D.verbose > 1
          fprintf('BrainAGE [years]:\n');
          fprintf('%20s\t','Group');
          for o = 1:n_groups
            fprintf('%20s\t',deblank(D.name_groups(o,:)));
          end
        
          fprintf('\n'); fprintf('%20s\t','Mean');   for o = 1:n_groups, fprintf('%20.3f\t',avg_BrainAGE(o)); end
          fprintf('\n'); fprintf('%20s\t','Median'); for o = 1:n_groups, fprintf('%20.3f\t',median_BrainAGE(o)); end
          fprintf('\n'); fprintf('%20s\t','SD');     for o = 1:n_groups, fprintf('%20.3f\t',SD_BrainAGE(o)); end
          fprintf('\n');
        end
  
        % ANOVA + T-Test
        if n_groups > 1
          if exist('BA_anova1')
            group = ones(length(D.ind_groups{1}),1);
            for o = 2:n_groups
              group = [group; o*ones(length(D.ind_groups{o}),1)];
            end
            if D.verbose > 1
              fprintf('\nANOVA P-value: ');
              for r = 1:D.n_regions
                Panova = BA_anova1(data(:,r),group,'off');
                fprintf('%f ',Panova);
              end
              fprintf('\n');
            end
            

            warning off
            for r = 1:D.n_regions
              P = zeros(n_groups,n_groups);
    
              for o = 1:n_groups
                for p = 1:n_groups
                  if o == p
                    P(o,p) = 0.5;
                  else
                    [~,P(o,p)] = BA_ttest2(BrainAGE(D.ind_groups{o},r),BrainAGE(D.ind_groups{p},r),0.05,'left');
                  end
                end
              end
                
              if D.n_regions > 1
                if is_surf
                  region_str = region_names_surf{r};
                else
                  region_str = region_names{r};
                end
              else
                region_str = '';
              end
              fprintf('T-test P-value (one-tailed) %s\n',region_str);
              fprintf('%20s\t','Group');
              for o = 1:n_groups
                fprintf('%20s\t',deblank(D.name_groups(o,:)));
              end
              fprintf('\n');
    
              for o = 1:n_groups
                fprintf('%20s\t',deblank(D.name_groups(o,:)));
                for p = 1:n_groups
                  if P(o,p) <= 0.001 || P(o,p) >= 0.999
                    fprintf('%20g***\t',P(o,p));
                  elseif P(o,p) <= 0.01 || P(o,p) >= 0.99
                    fprintf('%20g **\t',P(o,p));
                  elseif P(o,p) <= 0.05 || P(o,p) >= 0.95
                    fprintf('%20g  *\t',P(o,p));
                  else
                    fprintf('%20g\t',P(o,p));
                  end
                end
                fprintf('\n');
              end
          
              if any(P<=0.05)
                fprintf('****************************\n');
                fprintf('Significant result found\n');
                fprintf('****************************\n\n');
              end
            end
            warning on
                
          else
            fprintf('Warning: BA_anova1 not found.\n');
          end

 
        elseif isfield(D,'corr') % estimate correlation for one group
          for r = 1:D.n_regions
            [R, P] = corrcoef(BrainAGE(~isnan(D.corr),r),D.corr(~isnan(D.corr)));
            fprintf('Correlation r=%g (P=%g)\n',R(1,2),P(1,2));
          end
        end
        
      end
    end
  end        
end

% estimate ensembles
if multiple_BA

  if D.n_regions > 1
    BA_unsorted0 = reshape(BA_unsorted,D.n_data,D.n_regions,size(BA_unsorted,2)/D.n_regions);
    BA_unsorted = zeros(D.n_data,size(BA_unsorted,2)/D.n_regions,D.n_regions);
    for k = 1:D.n_data
      BA_unsorted(k,:,:) = squeeze(BA_unsorted0(k,:,:))';
    end
  end

  if isfield(D,'ensemble_method') && D.ensemble_method == 4, fprintf('GPR ensemble approach is only useful for k-fold validation.\n'); end
  BA_unsorted_weighted  = ensemble_models(BA_unsorted,D.age_test,D);
    
  if D.verbose > 0 && D.trend_degree >= 0
    fprintf('\n===========================================================\n'); 
    str_trend = {'No age correction','Age correction using BA','Age correction using PredicatedAge (Cole)'};
    co = 0:2;
    co(co == D.trend_method) = [];
    for i=co
      D1 = D;
      D1.trend_method = i;
      BA_unsorted_weighted1 = apply_trend_correction(BA_unsorted_weighted,age,D1);
      fprintf('\n%s:\n',str_trend{i+1});

      MAE_weighted = mean(abs(BA_unsorted_weighted1(D.ind_adjust,:)));
      cc = corrcoef([BA_unsorted_weighted1+age age]);
      fprintf('Weighted MAE/correlation (testdata): '); 
      for r = 1:D.n_regions
        fprintf(' %2.3f/%1.3f',MAE_weighted(r),cc(end,r));
      end
      fprintf('\n');
    end
    fprintf('\n===========================================================\n'); 
  end

  % estimate MAE and correlation using BA without normalization
  BA_unsorted_weighted0 = BA_unsorted_weighted;
  % apply final trend correction to weighted model
  if D.trend_degree >= 0
    BA_unsorted_weighted0 = apply_trend_correction(BA_unsorted_weighted0,age,D);
  end
  
  MAE_ctl_weighted = mean(abs(BA_unsorted_weighted0(D.ind_adjust,:)));

  cc = corrcoef([BA_unsorted_weighted0+age age]);
  fprintf('-----------------------------------------------------------------\n'); 
  fprintf('Weighted MAE/correlation (testdata): '); 
  for r = 1:D.n_regions
    fprintf(' %2.3f/%1.3f',MAE_ctl_weighted(r),cc(end,r));
  end
  fprintf('\n');
  fprintf('-----------------------------------------------------------------\n'); 

  
  % apply BA normalization w.r.t. MAE to make BA less dependent from training
  % sample (size) and scale it to MAE of 5
  if D.normalize_BA
    fprintf('\nApply normalization to BrainAGE using MAE.\n');
    BA_unsorted_weighted1 = apply_trend_correction(BA_unsorted_weighted,age,D,0);
    BA_unsorted_weighted = D.normalize_BA * BA_unsorted_weighted ./ (mean(abs(BA_unsorted_weighted1(D.ind_adjust,:))));
  end

  % check correlation coefficient of unsorrected data that should be not tto low
  cc = corrcoef([BA_unsorted_weighted(D.ind_adjust)+age(D.ind_adjust) age(D.ind_adjust)]);
  if cc(1,2) < 0.4
    spm('alert*','Correlation between predicted and chronological age (without age correction) is too low to be meaningful.');
    fprintf('\n********************************************************************\n'); 
    warning('Correlation between predicted and chronological age (without age correction) is too low to be meaningful: cc = %1.3f!',cc(1,2));
    fprintf('********************************************************************\n\n'); 
  end  
  
  if size(BA_unsorted_weighted,2) == 1
    figure(22)
    clf
    hold on 
    [p,~] = polyfit(D.age_test(D.ind_adjust),D.age_test(D.ind_adjust)+BA_unsorted_weighted(D.ind_adjust),1);
    yfit = polyval(p,[0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))]);
  
    for i = 1:n_groups
      plot(D.age_test(D.ind_groups{i}),D.age_test(D.ind_groups{i})+BA_unsorted_weighted(D.ind_groups{i}),'*','color',groupcolor(i,:))
    end
    line([0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],[0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],...
      'Color',[0 0 0],'lineWidth',2);
    line([0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],yfit,...
      'Color',[1 0 0],'lineWidth',2);    hold off
    title('Data')
    xlabel('Chronological Age [years]')
    ylabel('Predicted Age [years]')
    legend(char(D.name_groups, 'Reference', 'Linear Age Fit'),'Location','NorthWest')
    set(gca,'FontSize',20);
    set(gcf,'MenuBar','none');
  end

  % apply final trend correction to weighted model
  if D.trend_degree >= 0
    BA_unsorted_weighted = apply_trend_correction(BA_unsorted_weighted,age,D);
  end

  fprintf('\n');
  BA_weighted = BA_unsorted_weighted(ind_groups,:);
  
  if D.verbose >1
    figure(23)
    plot(D.age_test(D.ind_adjust),D.age_test(D.ind_adjust)+BA_unsorted_weighted(D.ind_adjust),'.')
    hold on 
    line([0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],[0.9*min(D.age_test(D.ind_adjust)) 1.1*max(D.age_test(D.ind_adjust))],...
      'Color',[0 0 0],'lineWidth',2);
    hold off
    title('Test Data')
    xlabel('Chronological Age [years]')
    ylabel('Predicted Age [years]')
    set(gca,'FontSize',20);

    figure(24)
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
    set(gcf,'MenuBar','none');

  end
  
  if n_groups > 1
    if exist('BA_anova1')

      if D.verbose > 1
        fprintf('\nWeighted ANOVA P-value: ');
        for r = 1:D.n_regions
          Panova = BA_anova1(BA_weighted(:,r),group,'off');
          fprintf('%f ',Panova);
        end
        fprintf('\n');
      end
    
      for r = 1:D.n_regions
        P = zeros(n_groups,n_groups);
        for o = 1:n_groups
          for p = 1:n_groups
            [~,P(o,p)] = BA_ttest2(BA_unsorted_weighted(D.ind_groups{o},r),BA_unsorted_weighted(D.ind_groups{p},r),0.05,'left');
          end
        end
      
        if D.n_regions > 1
          if is_surf
            region_str = region_names_surf{r};
          else
            region_str = region_names{r};
          end
        else
          region_str = '';
        end
        fprintf('T-test P-value (one-tailed) %s\n',region_str);
        fprintf('%20s\t','Group');
        for o = 1:n_groups
          fprintf('%20s\t',deblank(D.name_groups(o,:)));
        end
        fprintf('\n');
      
        for o = 1:n_groups
          fprintf('%20s\t',deblank(D.name_groups(o,:)));
          for p = 1:n_groups
            if P(o,p) <= 0.001 || P(o,p) >= 0.999
              fprintf('%20g***\t',P(o,p));
            elseif P(o,p) <= 0.01 || P(o,p) >= 0.99
              fprintf('%20g **\t',P(o,p));
            elseif P(o,p) <= 0.05 || P(o,p) >= 0.95
              fprintf('%20g  *\t',P(o,p));
            else
              fprintf('%20g\t',P(o,p));
            end
          end
          fprintf('\n');
        end
      
        if ~isempty(find(P<=0.05, 1))
          fprintf('****************************\n');
          fprintf('Significant result found\n');
          fprintf('****************************\n\n');
        end
      end
    else
      fprintf('BA_anova1 not found.\n');
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
if multiple_BA
      
  ind_groups = [];
  for o = 1:n_groups
    ind_groups = [ind_groups; D.ind_groups{o}];
    data_cell{o} = BA_unsorted_weighted(D.ind_groups{o},:);
    avg_BrainAGE(o,:) = mean(BA_unsorted_weighted(D.ind_groups{o},:));
    age_test(o,1) = mean(D.age_test(D.ind_groups{o}));
    age_test(o,2) = std(D.age_test(D.ind_groups{o}));
    median_BrainAGE(o,:) = median(BA_unsorted_weighted(D.ind_groups{o},:));
    SD_BrainAGE(o,:) = std(BA_unsorted_weighted(D.ind_groups{o},:));
  end
    
  % print BrainAGE of groups
  fprintf('%20s\t','Group');
  for o = 1:n_groups
    fprintf('%20s\t',deblank(D.name_groups(o,:)));
  end

  if D.n_regions > 1, fprintf('\n'); end
  for r = 1:D.n_regions
    if D.n_regions > 1
      if is_surf
        fprintf('Region: %s\n',region_names_surf{r});
      else
        fprintf('Region: %s\n',region_names{r});
      end
    else
      fprintf('\n'); 
    end
    fprintf('%20s\t','Mean');   for o = 1:n_groups, fprintf('%20.3f\t',avg_BrainAGE(o,r)); end
    fprintf('\n'); fprintf('%20s\t','Median'); for o = 1:n_groups, fprintf('%20.3f\t',median_BrainAGE(o,r)); end
    fprintf('\n'); fprintf('%20s\t','SD');     for o = 1:n_groups, fprintf('%20.3f\t',SD_BrainAGE(o,r)); end
    fprintf('\n');
  end

  if D.verbose && sum(isnan(BA_unsorted_weighted(:))) == 0
    f = figure(24);
    set(f, 'Position',[10 10 900 800])
    
    % show spiderplot for regional BrainAGE
    if D.parcellation
        combined_BA = zeros(numel(data_cell),size(BA_unsorted_weighted,2));
        for l = 1:numel(data_cell)
            combined_BA(l,:) = feval(D.spiderplot.func,BA_unsorted_weighted(D.ind_groups{l},:));
        end
        if isfield(D,'spiderplot')
          if isfield(D.spiderplot,'range')
          BA_spider_plot(combined_BA, 'Names', D.name_groups, 'Colors', groupcolor, 'Parent', f, 'range', D.spiderplot.range);
          else
            BA_spider_plot(combined_BA, 'Names', D.name_groups, 'Colors', groupcolor, 'Parent', f);
          end
        end
        if strcmpi(D.spiderplot.func,'median')
          set(f,'Name','Median Weighted BrainAGE','MenuBar','none');
        elseif strcmpi(D.spiderplot.func,'mean')
          set(f,'Name','Mean Weighted BrainAGE','MenuBar','none');
        end
    else % otherwise use boxplot
        cat_plot_boxplot(data_cell,opt);
        if style == 1
          set(gca,'XTick',1:n_groups,'XTickLabel',D.name_groups);
          ylabel('BrainAGE [years]');
        else
          set(gca,'YTick',1:n_groups,'YTickLabel',D.name_groups(n_groups:-1:1,:));
          xlabel('BrainAGE [years]');
        end
        set(f,'Name','Weighted BrainAGE','MenuBar','none');
    end

    set(gca,'FontSize',20);
  end
end

if (D.age_range(1)-min(D.age_test)) > 2 || (max(D.age_test)-D.age_range(2)) > 2
  fprintf('\n********************************************************************\n'); 
  fprintf('Warning: Defined age range of training sample (%3.1f..%3.1f years) differs from real age range of sample %g..%g\n',...
      D.age_range(1),D.age_range(2),min(D.age_test),max(D.age_test));
  fprintf('********************************************************************\n'); 
end


%-------------------------------------------------------------------------------
function [BA_weighted, PredictedAge_weighted] = ensemble_models(BA, age_all, D, ind_test_array, ind_train_array)
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
  BA_weighted = squeeze(BA);
  if isfield(D,'define_cov')
    PredictedAge_weighted = BA_weighted;
  else
    PredictedAge_weighted = BA_weighted + age;
  end
  fprintf('Disable model ensemble because we just have one model.\n');
  return
end

% use weighted mean as default if no ensemble method is given
if isfield(D,'ensemble')
  ensemble_method = D.ensemble;
else
  ensemble_method = 5;
end

if (nargin < 4) && (ensemble_method == 4)
  fprintf('For GPR ensemble method you have to define 5 arguments.\n');
  ensemble_method = 1;
end

% remove columns with NaNs where model estimation failed
BA_corrected = BA;
n_regions = size(BA,3);
indy_all = [];
for r = 1:n_regions
  [~,indy] = find(isnan(BA_corrected(:,:,r)));
  if ~isempty(indy)
    indy_all = [indy_all unique(indy)];
  end
end
if ~isempty(indy_all)
  BA_corrected(:,unique(indy_all),:) = [];
end

if isfield(D,'define_cov')
  PredictedAge_corrected = BA_corrected;
else
  PredictedAge_corrected = BA_corrected + age;
end

switch ensemble_method
case 0   % use model with lowest MAE
  
  BA_weighted = zeros(size(PredictedAge_corrected,1),D.n_regions);
  for r = 1:n_regions
    BA_ind = BA_corrected(D.ind_adjust,:,r);
    [~, mn_ind] = min(mean(abs(BA_ind)));
  
    fprintf('\nModel with lowest MAE for region %d is %d', r, mn_ind);
  
    BA_weighted(:,r) = PredictedAge_corrected(:,mn_ind,r) - age;
  end
  fprintf('\n');

case 1   % use GLM estimation to minimize MAE
  
  BA_weighted = zeros(size(PredictedAge_corrected,1),D.n_regions);
  fprintf('Estimated weights using %d subjects:\n',numel(D.ind_adjust));
  for r = 1:D.n_regions
    PredictedAge_ind = PredictedAge_corrected(D.ind_adjust,:,r);
    Beta = pinv(PredictedAge_ind)*age(D.ind_adjust);
  
    fprintf('%.2f ',Beta);
    fprintf('\n');
  
    weighted_PredictedAge = PredictedAge_corrected(:,:,r)*Beta;
    BA_weighted(:,r) = weighted_PredictedAge - age;
  end

case 2 % simple mean of all models
  
  BA_weighted = squeeze(mean(BA_corrected,2));
case {3, -3}   % use GLM estimation to maximize group differences or correlation (EXPERIMENTAL!)

  if ~isfield(D,'contrast')
    error('D.contrast has to be defined.');
  end

  % F-contrast ?
  if size(D.contrast,1) == numel(D.ind_groups) && size(D.contrast,2) == numel(D.ind_groups)
    c = D.contrast;
    D.contrast = [];
    ind = [];
    for i = 1:numel(D.ind_groups)
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
  
  disp('#############################################');
  warning('Experimental: Not yet working properly!');
  disp('#############################################');

  if group_diff
    X = []; Y = [];
    for i = 1:size(D.contrast,2)
      if any(D.contrast(:,i) ~= 0)
        X = [X; D.contrast(:,i)*ones(numel(D.ind_groups{i}),1)];  
        Y = [Y; mean(BA_corrected(D.ind_groups{i},:,:),3)];  
      end
    end
  else
    X = D.contrast;
    Y = mean(BA_corrected,3);
  end
  
  % use only those data that are defined in ind_groups
  ind = [];
  for i = 1:numel(D.ind_groups)
    ind = [ind; D.ind_groups{i}];
  end
    
  if group_diff
    Yind = Y;
    Xind = X;
  else
    % we have to excluded NaNs for Beta estimation
    ind_finite = all(isfinite(X),2);
    Yind = Y(ind(ind_finite),:);
    Xind = X(ind_finite,:);
  end
  
  if ensemble_method > 0
    Beta = pinv(Yind)*Xind;
  else % use non-linear optimization
    Beta = nonlin_optim(Yind, Xind);
  end


  fprintf('\nPredicted weights using %d subjects:\t',numel(D.ind_adjust));
  fprintf('%.2f ',Beta);
  fprintf('\n');

  n_regions = size(BA_corrected,3);
  BA_weighted = zeros(size(BA_corrected,1), n_regions);
  for i=1:n_regions
    BA_weighted(:,i) = sum(BA_corrected(:,:,i)*Beta,2);
  end

  % scale estimated BA values by ratio between SD of original and estimated BA values to get the same range
  BA_weighted = BA_weighted*mean(std(BA_corrected(:)))/mean(std(BA_weighted));
case 4   % Stacking: use GPR to combine models
        
  PredictedAge_ind = PredictedAge_corrected;
  BA_weighted = zeros(size(PredictedAge_corrected,1),D.n_regions);
  
  for i = 1:numel(ind_test_array)
    
    % get indices for training and test and only include data that are
    % defined in D.ind_adjust
    ind_train = ind_train_array{i};
    ind_train = ind_train(ismember(ind_train,D.ind_adjust));
    ind_test  = ind_test_array{i};
    ind_test  = ind_test(ismember(ind_test,D.ind_adjust));
    
    for r = 1:D.n_regions
      Y_train = PredictedAge_ind(ind_train,:,r);
      Y_test  = PredictedAge_ind(ind_test,:,r);
    
      % multivariate regression using GPR
      PredictedAge = BA_gpr_core(Y_train, age(ind_train), Y_test, 10, 1);
      BA_weighted(ind_test,r) = PredictedAge-age(ind_test);
    end
  end
  
%  BA_weighted = BA_weighted';

case 5 % weighted average of all models w.r.t. MAE
  BA_weighted = zeros(size(PredictedAge_corrected,1),D.n_regions);
  for r = 1:D.n_regions
    MAE = mean(abs(BA_corrected(:,:,r)));
    weight = 1./MAE.^2;
    weight = weight/mean(weight);
    PredictedAge_weighted = mean(PredictedAge_corrected(:,:,r).*weight,2);
    BA_weighted(:,r) = PredictedAge_weighted - age;
  end

case {6, -6}   % use GLM estimation for mean tissue to maximize group differences or correlation (EXPERIMENTAL!)

  if ~isfield(D,'contrast')
    error('D.contrast has to be defined.');
  end

  % number of tissue types
  n_tissues = numel(D.seg_array);

  if n_tissues < 2
    error('D.seg_array must have at least two entries (i.e. rp1/rp2).');
  end
  
  % F-contrast ?
  if size(D.contrast,1) == numel(D.ind_groups) && size(D.contrast,2) == numel(D.ind_groups)
    c = D.contrast;
    D.contrast = [];
    ind = [];
    for i = 1:numel(D.ind_groups)
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
  
  disp('#############################################');
  warning('Experimental: Not yet working properly!');
  disp('#############################################');

  if group_diff
    X = []; Y = [];
    for i = 1:size(D.contrast,2)
      if any(D.contrast(:,i) ~= 0)
        X = [X; D.contrast(:,i)*ones(numel(D.ind_groups{i}),1)];
        % use mean of tissue types (every n_tissues entry) to build model with just n_tissues weights
        mean_tissue = [];
        for j = 1:n_tissues
          mean_tissue = [mean_tissue mean(squeeze(mean(BA_corrected(D.ind_groups{i},j:n_tissues:end,:),2)),2)];
        end
        Y = [Y; mean_tissue];  
      end
    end
  else
    X = D.contrast;
    mean_tissue = [];
    for j = 1:n_tissues
      mean_tissue = [mean_tissue mean(squeeze(mean(BA_corrected(:,j:n_tissues:end,:),2)),2)];
    end
    Y = mean_tissue;
  end
  
  mean_tissue = [];
  for j = 1:n_tissues
    mean_tissue = [mean_tissue mean(BA_corrected(:,j:n_tissues:end,:),2)];
  end
  Yall = mean_tissue;
  
  % use only those data that are defined in ind_groups
  ind = [];
  for i = 1:numel(D.ind_groups)
    ind = [ind; D.ind_groups{i}];
  end
  
  if group_diff
    Yind = Y;
    Xind = X;
  else
    % we have to excluded NaNs for Beta estimation
    ind_finite = all(isfinite(X),2);
    Yind = Y(ind(ind_finite),:);
    Xind = X(ind_finite,:);
  end
  
  if ensemble_method > 0
    Beta = pinv(Yind)*Xind;
  else % use non-linear optimization
    Beta = nonlin_optim(Yind, Xind);
  end

  fprintf('\nPredicted weights using %d subjects:\t',numel(D.ind_adjust));
  fprintf('%.2f ',Beta);
  fprintf('\n');

  n_regions = size(BA_corrected,3);
  BA_weighted = zeros(size(BA_corrected,1), n_regions);
  for i=1:n_regions
    BA_weighted(:,i) = sum(Yall(:,:,i)*Beta,2);
  end

  % scale estimated BA values by ratio between SD of original and estimated BA values to get the same range
  BA_weighted = BA_weighted*mean(std(BA_corrected(:)))/mean(std(BA_weighted));
  
end

if isfield(D,'define_cov')
  PredictedAge_weighted = BA_weighted;
else
  PredictedAge_weighted = BA_weighted + age;
end

%-------------------------------------------------------------------------------
function [BrainAGE, PredictedAge, Adjustment] = apply_trend_correction(BrainAGE, age, D, verbose)
%-------------------------------------------------------------------------------
% BrainAGE     - adjusted BrainAGE
% PredictedAge - adjusted PredictedAge
% Adjustment   - array of adjustments for indexed data

PredictedAge = BrainAGE + age;
if D.trend_degree < 0 || ~D.trend_method
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
  for j = 1:length(D.ind_adjust)
    if site_adjust(D.ind_adjust(j)) == i
      ind_site_adjust = [ind_site_adjust; D.ind_adjust(j)];
    end
  end

  % apply trend correction for indexed data
  if ~isempty(ind_site_adjust)
  
    % test whether sample size is too small
    if length(ind_site_adjust) < 15
      warning('Sample #%d for site-specific trend-correction is too small (n=%d). Use only offset-correction by Mean instead.\n',...
          i,length(ind_site_adjust));
      offset = median(BrainAGE(ind_site_adjust,:));
      BrainAGE(ind_site,:) = BrainAGE(ind_site,:) - offset;
      PredictedAge(ind_site,:) = PredictedAge(ind_site,:) - offset;
      
    % test whether age range is too small
    elseif (max(age(ind_site_adjust)) - min(age(ind_site_adjust))) < 5
      warning('Age range for site-specific trend-correction is too small (%g). Use only offset-correction by Mean instead.\n',...
          (max(age(ind_site_adjust)) - min(age(ind_site_adjust))));
      offset = median(BrainAGE(ind_site_adjust,:));
      BrainAGE(ind_site,:) = BrainAGE(ind_site,:) - offset;
      PredictedAge(ind_site,:) = PredictedAge(ind_site,:) - offset;
      
    % everything else
    else

      % use predicted age for obtaining trend
      if D.trend_method == 2
        Y = PredictedAge;
      else
        % use BrainAGE for obtaining trend
        Y = BrainAGE;
      end
            
      % polynomial trend
      G = ones(length(age(ind_site)),1);
      G_indexed = ones(length(age(ind_site_adjust)),1);

      % polynomial expansion
      if D.trend_degree > 0
        G = [G cat_stat_polynomial(age(ind_site),D.trend_degree)];
        G_indexed = [G_indexed cat_stat_polynomial(age(ind_site_adjust),D.trend_degree)];
      end

      if verbose, fprintf('Remove trend degree %d using %d subjects of site %d.\n',D.trend_degree,length(ind_site_adjust),i); end
        
      % estimate beta only for indexed data (e.g. control subjects)
      Beta = pinv(G_indexed)*Y(ind_site_adjust,:);
      Beta = mean(Beta,2);

      % and remove effects for all data
      GBeta = G*Beta;

      BrainAGE0 = BrainAGE;
      
      % use predicted age for obtaining trend
      if D.trend_method == 2
        PredictedAge(ind_site,:)  = (PredictedAge(ind_site,:) - Beta(1));
        if D.trend_degree == 1
          PredictedAge(ind_site,:) = PredictedAge(ind_site,:)/Beta(2);
        end
      else
        % use BrainAGE for obtaining trend
        PredictedAge(ind_site,:)  = PredictedAge(ind_site,:)  - GBeta;
      end

      BrainAGE = PredictedAge - age;
      Adjustment{i} = BrainAGE0(ind_site,:) - BrainAGE(ind_site,:);
    end
  else
    if D.trend_degree >= 0
      if verbose, fprintf('Warning: No subjects found in sample #%d for site-specific trend-correction\n',i); end
    end
  end
end

avg_BrainAGE = mean(BrainAGE(D.ind_adjust,:,:), 1);
BrainAGE = BrainAGE - avg_BrainAGE;
PredictedAge = PredictedAge - avg_BrainAGE;

MAE = mean(abs(BrainAGE));

if MAE0/MAE > 4
  fprintf('Warning: Large discrepancy between MAE before and after correction which points to a too narrow age range of training data!\n');
end

%-------------------------------------------------------------------------------
function C = nejm
%-------------------------------------------------------------------------------
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
C = reshape(sscanf(C(:,2:end)','%2x'),3,[])./255;
C = C';

function Beta = nonlin_optim(Y, X)
% Use non-linear optimization to solve the equivalent of pinv(Y)*X, but with the
% constrain that weights (betas) should be positive and in the range 0.001..0.999
% Requirement: Optimization Toolbox

num_ensembles = size(Y,2);

% Objective function: Mean Squared Error of weighted predictions
objective = @(weights) mean((X - Y * weights).^2);

% Initial guess
initial_weights = ones(num_ensembles, 1) / num_ensembles;

% Linear equality constraint to make weights sum to 1
Aeq = ones(1, num_ensembles);
beq = 1;

% Bounds to ensure weights are non-negative
lb = zeros(num_ensembles, 1) + 0.001;
ub = ones(num_ensembles, 1) - 0.001;

% Options: Display iterations
options = optimoptions('fmincon', 'Display', 'iter');

% Optimization
Beta = fmincon(objective, initial_weights, [], [], Aeq, beq, lb, ub, [], options);

