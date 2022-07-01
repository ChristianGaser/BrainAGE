function [BA_all, P_all, BA_unsorted_all, MAE_all, Names_all, BA_weighted, BA_unsorted_weighted, MAE_weighted] = cg_BrainAGE_ui(D)
% [BA_all, P_all, BA_unsorted_all, MAE_all, Names_all, BA_weighted, BA_unsorted_weighted, MAE_weighted] = cg_BrainAGE_ui(D)
% user interface for BrainAGE estimation
%
% D.training_sample - training sample
%     Healthy adults:
%         IXI547    - IXI-database
%         IXI410    - IXI-database (1:4:547 removed)
%         OASIS316  - OASIS
%         Ship1000  - SHIP (internal use only)
%     ADNI231Normal - ADNI Normal-sc-1.5T 
%Long_Australia2611 - PATH study Australia (internal use only)
%
%     Children data:
%         NIH394    - NIH objective 1
%         NIH876    - NIH release 4.0
%         NIH755    - NIH release 5.0
%
%     Children + adults:
%         fCONN772 - fcon-1000 (8-85 years)
%
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
%                     [0 Inf] use all data
%                     [50 80] use age range of 50..80
%                     if not defined use min/max of age of test data
% D.nuisance        - additionally define nuisance parameter for covarying out (e.g. gender)
% D.ind_adjust      - define indices for adjusting data according to trend defined with D.trend_degree
%                     usually this is the control group that is used in order to adjust the data
% D.site_adjust     - if data are acquired at different sites (e.g. using different scanners or sequences) the trend 
%                     correction should be estimated and applied for each site seperately
% D.ind_groups      - define indices of groups, if D.ind_adjust is not given, then the first group of this index will
%                     be used for adjusting data according to trend defined with D.trend_degree
% D.ind_train       - define indices of subjects used for training (e.g. to limit training on male subjects)
% D.trend_degree    - estimate trend with defined order using healthy controls and apply it to all data (set to -1 for skipping trend correction)
% D.n_fold          - n-fold validation if training and test sample are the same or only one is defined (10-fold as default)
% D.validate        - Apply n-fold validation for 'inner_loop' (as default) or 'outer_loop'
%                     'inner_loop' means here that the weightings of the different models are estimated in nested loops inside the n-fold
%                     'outer_loop' estimates the weightings after n-fold validation is finished
% D.weighting       - weighting of different models
%                     1 use GLM estimation to minimize MAE
%                     2 simple median (as default)
%                     3 use GLM estimation to maximize group differences (EXPERIMENTAL!)
% D.contrast        - define contrast to maximize group differences (use only if D.weighting=3) (e.g. [0 1 -1])
% D.dir             - directory for databases and code
% D.verbose         - verbose level (default=1), set to "0" to suppress long outputs
% D.threshold_std   - all data with a standard deviation > D.threshold_std of mean covariance are excluded
%                     (after covarying out effects of age)
% D.corr            - additionally define parameter that can be correlated to BrainAGE if only one group is given
% D.style           - plot-style: 1: old style with vertical violin-plot; 2: new style with horizontal density plot
% D.groupcolor      - matrix with (group)-bar-color(s), use jet(numel(data)) or other color functions
%
% Parameter search
% ---------------
% Some selected parameters can be also defined as ranges to try different parameter settings.
% Examples:
% D.trend_degree = [0 1 2];
% D.threshold_std = [Inf 1 2];
% D.age_range = {[20 50],[20 60],[20 70]};
% D.res_array    = {'4','8'};   
% D.smooth_array = {'s4','s8'}; 
% or
% age = 50:10:70;
% for i=1:numel(age)
%   D.age_range{i} = [20 age(i)];
% end
% D.weighting = 1; % minimize MAE
%
% Output
% -------
% BA_all            - cell array of BrainAGE values sorted by group definitions in D.ind_groups 
% P_all             - cell array of P-values of One-Way ANOVA
% BA_unsorted_all   - cell array of unsorted (originally ordered) BrainAGE values
% MAE_all           - cell array of mean absolute errors (additionally for each group)
% Names_all         - cell array of names
% weighted outputs:
% BA_weighted           - weighted BrainAGE sorted by group definitions in D.ind_groups
% BA_unsorted_weighted  - weighted BrainAGE
% MAE_weighted          - MAE of weighted BrainAGE
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_BrainAGE_ui.m 4 2015-08-28 08:57:46Z gaser $

if ~isfield(D,'trend_degree')
  trend_degree = 2;
else
  trend_degree = D.trend_degree;
end

if ~isfield(D,'style')
  style = 2;
else
  style = D.style;
end

if ~isfield(D,'threshold_std')
  threshold_std = Inf;
else
  threshold_std = D.threshold_std;
end

% convert all to cell format to allow multiple entries for lower/upper age ranges
if isfield(D,'age_range')
  if iscell(D.age_range)
    age_range = D.age_range;
  else
    age_range{1} = D.age_range;
  end
  n_age_range = numel(age_range);
else
  n_age_range = 1;
end

% this is just for compatbility with older scripts
if isfield(D,'seg') & ~isfield(D,'seg_array')
  D.seg_array = D.seg;
  D = rmfield(D,'seg');
end

% convert to cell if necessary
if ~iscell(D.seg_array)
  D.seg_array = cellstr(D.seg_array);
end

if ~iscell(D.smooth_array)
  D.smooth_array = cellstr(D.smooth_array);
end

if ~iscell(D.res_array)
  D.res_array = cellstr(D.res_array);
end

% verbose level
if ~isfield(D,'verbose')
  D.verbose = 1;
end

% fill the missing field
if ~isfield(D,'data')
  D.data = D.training_sample{1};
end

% load first data set to get data size
if ~isfield(D,'n_data')
  ind_plus = strfind(D.seg_array{1},'+');
  if ~isempty(ind_plus)
    seg_array = D.seg_array{1}(1:ind_plus-1);
  else
    seg_array = D.seg_array{1};
  end
  name = [D.smooth_array{1} seg_array '_' D.res_array{1} 'mm_' D.data D.relnumber];
  load(name);
  D.n_data = size(Y,1);
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

BA_unsorted_weighted = []; MAE_weighted = [];
BA_all = []; P_all = []; BA_unsorted_all = [];
Names_all = []; BA_weighted = []; MAE_all = [];

% run n-fold validation if no data field is given or validation with n_fold is defined
if ((~isfield(D,'data') | ~isfield(D,'training_sample')) | isfield(D,'n_fold')) & ~isfield(D,'run_validation')
  % use inner_loop n-fold validation as default
  if ~isfield(D,'validate')
    D.validate = 'inner_loop';
  end
    
  % use 10-fold as default
  if ~isfield(D,'n_fold')
    D.n_fold = 10;
  end
  
  % indicate that validation is running and will not called in nested loops
  D.run_validation = 1;
  
  EA_all   = zeros(D.n_data,1);
  ind_all  = [];
  age_all  = [];
  
  for j=1:D.n_fold
    ind = j:D.n_fold:D.n_data;
    D.ind_groups  = {ind};
    age_fold = age(ind);
    age_all = [age_all; age_fold];
    ind_all  = [ind_all ind];
    D.ind_adjust = ind;
    
    % build training sample using remaining subjects
    ind_train = 1:D.n_data;
    ind_train(ind) = [];
    D.ind_train = ind_train;

    [data, ~, ~, ~, Names_all] = cg_BrainAGE_ui(D);
    if j==1
      data_all = data;
    else
      for i=1:numel(data)
        data_all{i} = [data_all{i}; data{i}];
      end
    end
  
    if strcmp(D.validate,'inner_loop')    
      EA_all(ind) = get_EstimatedAge(data,age_fold);
    end
    
  end
  
  if strcmp(D.validate,'inner_loop')
    BA_unsorted_weighted = EA_all-age;
  else
    EstimatedAge_unsorted_weighted = get_EstimatedAge(data_all,age_all);
    BA_unsorted_weighted = EstimatedAge_unsorted_weighted-age_all;
    BA_unsorted_weighted(ind_all) = BA_unsorted_weighted;
  end

  % apply trend correction
  BA_unsorted_weighted = apply_trend_correction(BA_unsorted_weighted,age,D);

  MAE_weighted = mean(abs(BA_unsorted_weighted));
  fprintf('weighted MAE = %g\n\n',MAE_weighted);
  
  return
end

% check whether additional fields for weighted BA are necessary
array_size = numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array)*n_age_range*numel(threshold_std)*numel(trend_degree);
multiple_BA = (array_size > 1);

% prepare output
BA_all          = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array),n_age_range,numel(threshold_std),numel(trend_degree));
BA_unsorted_all = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array),n_age_range,numel(threshold_std),numel(trend_degree));
P_all           = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array),n_age_range,numel(threshold_std),numel(trend_degree));
Names_all       = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array),n_age_range,numel(threshold_std),numel(trend_degree));
EA_unsorted     = [];
EA              = [];

count  = 0;
count0 = 0;

for i = 1:numel(D.res_array)
  for j = 1:numel(D.smooth_array)
    for o = 1:numel(D.seg_array)
  
      count = count + 1;
    
      D.res = D.res_array{i};
      D.smooth = D.smooth_array{j};
      seg = D.seg_array{o};
      
      % remove old D.seg field if exist
      if isfield(D,'seg')
        D = rmfield(D,'seg');
      end
      
      % find potenial "+" indicating to combine segmentations
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
      name = [D.smooth D.seg{1} '_' D.res 'mm_' D.data D.relnumber];
      load(name);
  
      if D.verbose, fprintf('\n-----------------------------------------------------------------\n%s\n',name); end
      D.Y_test = Y;
      D.age_test = age;
  
      if ~isfield(D,'age_range')
        age_range{1} = [min(D.age_test) max(D.age_test)];
      end
      
      % use additional segmentation if defined
      for l = 2:n_seg
        name = [D.smooth D.seg{l} '_' D.res 'mm_' D.data D.relnumber];
        if D.verbose, fprintf('%s\n',name); end
        load(name);
        D.Y_test = [D.Y_test Y];
      end
      
      if D.verbose, fprintf('\n'); end
      if ~isfield(D,'ind_groups')
        D.ind_groups = {1:length(D.age_test)};
      end
      
      n_groups = numel(D.ind_groups);
      
      if ~isfield(D,'groupcolor')
        groupcolor = flipud(hsv(max(5,n_groups)));
      else
        groupcolor = D.groupcolor;
      end
  
      if count == 1
        MAE_all    = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array),n_age_range,numel(threshold_std),numel(trend_degree),n_groups);
        MAE_ctl    = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array),n_age_range,numel(threshold_std),numel(trend_degree));
      end
      
      if ~isfield(D,'nuisance')
        D.nuisance = [];
      end
      
      % build index for test data
      ind_test = [];
      for k = 1:n_groups
        if size(D.ind_groups{k},1) < size(D.ind_groups{k},2)
          D.ind_groups{k} = D.ind_groups{k}';
        end
        ind_test = [ind_test; D.ind_groups{k}];
      end
      
      for l=1:numel(threshold_std)
        for m=1:n_age_range
          for n=1:numel(trend_degree)
              
            count0 = count0 + 1;
            D.threshold_std = threshold_std(l);
            D.trend_degree = trend_degree(n);
            
            if ~isfield(D,'age_range')
              D.age_range = [min(D.age_test(ind_test)) max(D.age_test(ind_test))];
            else
              D.age_range = age_range{m};
              if numel(D.age_range) ~=2
                error('Age range has to be defined by two values (min/max)');
              end
            end
  
            % give warning if age range does not fit
            if D.age_range(1)>min(ceil(D.age_test(ind_test))) | D.age_range(2)<max(floor(D.age_test(ind_test)))
              fprintf('Defined age range of training sample min=%g..max=%g does not fit to real age range of sample %g..%g\n',...
                  D.age_range(1),D.age_range(2),min(D.age_test(ind_test)),max(D.age_test(ind_test)));
            end
            
            calc_LI = 0;
        
            % estimate lateralization index LI
            if isfield(D,'hemi')
              if strcmp(D.hemi,'li')
                calc_LI = 1;
              end
            end
      
            % estimate lateralization index LI
            if calc_LI
              D.hemi = 'lh';
              BrainAGE_lh = cg_BrainAGE(D);
              D.hemi = 'rh';
              BrainAGE_rh = cg_BrainAGE(D);
              BrainAGE = (BrainAGE_lh - BrainAGE_rh)./(BrainAGE_lh + BrainAGE_rh);
            else
              BrainAGE = cg_BrainAGE(D);
            end
            
            BrainAGE = apply_trend_correction(BrainAGE,D.age_test,D);
                      
            data_cell = cell(1,n_groups);
            data = [];
            avg_BrainAGE = zeros(1,n_groups);
            age_test = zeros(n_groups,2);
            median_BrainAGE = zeros(1,n_groups);
        
            MAE_ctl{count,m,l,n} = mean(abs(BrainAGE(D.ind_adjust)));
            ind_groups = [];
            allow_violin = 1;
            for k = 1:n_groups
              if length(D.ind_groups{k}) < 10
                allow_violin = 0;
              end
              ind_groups = [ind_groups; D.ind_groups{k}];
              data_cell{k} = BrainAGE(D.ind_groups{k});
              avg_BrainAGE(k) = mean(BrainAGE(D.ind_groups{k}));
              age_test(k,1) = mean(D.age_test(D.ind_groups{k}));
              age_test(k,2) = std(D.age_test(D.ind_groups{k}));
              median_BrainAGE(k) = median(BrainAGE(D.ind_groups{k}));
              data = [data; BrainAGE(D.ind_groups{k})];
              MAE_all{count,m,l,n,k} = mean(abs(BrainAGE(D.ind_groups{k})));
            end
        
            % save data for all spatial resolutions and smoothing sizes
            BA_unsorted_all{count0} = BrainAGE;
            BA_all{count,m,l,n}     = data;
  
            EA_unsorted = [EA_unsorted, BrainAGE + D.age_test];
            EA          = [EA, BrainAGE(ind_groups) + D.age_test(ind_groups)];
            
            Names_all{count,m,l,n} = [seg '-' D.smooth '-' D.res 'mm-std' num2str(D.threshold_std) '-trend' num2str(D.trend_degree)...
                 '-age' num2str(D.age_range(1)) '-' num2str(D.age_range(2))];
  
            if style == 1
              opt = struct('violin',2,'showdata',1,'groupcolor',groupcolor);
            else
              opt = struct('style',3,'groupcolor',groupcolor);
            end
            
            if ~allow_violin
              opt = struct('style',0,'groupcolor',groupcolor,'showdata',1);
              style = 1;
            end
            
            if D.verbose
              figure(21)
              tmp = cat_plot_boxplot(data_cell,opt);
        
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
              for k = 1:n_groups
                fprintf('%20s\t',deblank(D.name_groups(k,:)));
              end
            
              fprintf('\n'); fprintf('%20s\t','Mean'); for k = 1:n_groups, fprintf('%20.3f\t',age_test(k,1)); end
              fprintf('\n'); fprintf('%20s\t','SD');   for k = 1:n_groups, fprintf('%20.3f\t',age_test(k,2)); end
              fprintf('\n%20s\t','Sample size');
  
              for k = 1:n_groups
                fprintf('%20s\t',sprintf('n=%d',length(D.ind_groups{k})));
              end
              fprintf('\n\n');
  
            end
  
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
            % print BrainAGE of groups
            if D.verbose
              fprintf('%20s\t','Group');
              for k = 1:n_groups
                fprintf('%20s\t',deblank(D.name_groups(k,:)));
              end
            
              fprintf('\n'); fprintf('%20s\t','Mean');   for k = 1:n_groups, fprintf('%20.3f\t',avg_BrainAGE(k)); end
              fprintf('\n'); fprintf('%20s\t','Median'); for k = 1:n_groups, fprintf('%20.3f\t',median_BrainAGE(k)); end
              fprintf('\n');
        
            end 
      
            % ANOVA + T-Test
            if n_groups > 1
              if exist('cg_anova1')
                group = ones(length(D.ind_groups{1}),1);
                for k = 2:n_groups
                  group = [group; k*ones(length(D.ind_groups{k}),1)];
                end
                Panova = cg_anova1(data,group,'off');
                if D.verbose, fprintf('\nANOVA P-value: %f\n',Panova); end
      
                P = zeros(n_groups,n_groups);
      
                for o = 1:n_groups
                  for p = 1:n_groups
                    [H,P(o,p)] = cg_ttest2(BrainAGE(D.ind_groups{o}),BrainAGE(D.ind_groups{p}),0.05,'left');
                  end
                end
            
                P_all{count,m,l,n}  = Panova;
      
                fprintf('T-test P-value (one-tailed):\n');
                fprintf('%20s\t','Group');
                for k = 1:n_groups
                  fprintf('%20s\t',deblank(D.name_groups(k,:)));
                end
                fprintf('\n');
      
                for o = 1:n_groups
                  fprintf('%20s\t',deblank(D.name_groups(o,:)));
                  for p = 1:n_groups
                    if P(o,p) <= 0.001 | P(o,p) >= 0.999
                      fprintf('%20.7f***\t',P(o,p));
                    elseif P(o,p) <= 0.01 | P(o,p) >= 0.99
                      fprintf('%20.7f **\t',P(o,p));
                    elseif P(o,p) <= 0.05 | P(o,p) >= 0.95
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
                P_all = [];
              end
      
            elseif isfield(D,'corr') % estimate correlation for one group
              [R, P] = corrcoef(BA_all{count,m,l}(~isnan(D.corr)),D.corr(~isnan(D.corr)))
              P_all{count,m,l}  = P(1,2);     
            end
  
          end
        end
      end
    end
  end        
end

% estimate weightings
if multiple_BA

  % use averaging as default if no weighting method is given
  if isfield(D,'weighting')
    weighting_method = D.weighting;
  else
    weighting_method = 2;
  end
  
  switch weighting_method
  case 1   % use GLM estimation to minimize MAE
      
    EA_ind_adjust = EA_unsorted(D.ind_adjust,:);
    beta = EA_ind_adjust\D.age_test(D.ind_adjust);
    
    % estimate offset that remains after fitting
    offset = mean(D.age_test(D.ind_adjust)) - mean(EA_ind_adjust*beta);
    
    BA_unsorted_weighted = EA_unsorted*beta - D.age_test + offset;
    
    if 0 & D.verbose
      for i=1:numel(beta)
        fprintf('%s:\tweighting %g\n',Names_all{i},beta(i));
      end
    end
  case 2 % simple median of different models
    BA_unsorted_weighted = median(reshape(cell2mat(BA_unsorted_all),numel(BA_unsorted_all{1,1,1,1}), array_size),2);
  case 3   % use GLM estimation to maximize group differences (EXPERIMENTAL!)

    if ~isfield(D,'contrast')
      error('D.con has to be defined.');
    end

    if numel(D.contrast) ~= n_groups
      error('D.con has different size than number of groups.');    
    end

    warning('Experimental: Not yet working properly!');
    BA_unsorted = EA_unsorted - D.age_test;

    contrast = [0 1 -1];
    X = []; Y = [];
    for i=1:numel(D.contrast)
      if D.contrast(i) ~= 0
        X = [X; D.contrast(i)*ones(numel(D.ind_groups{i}),1)];  
        Y = [Y; BA_unsorted(D.ind_groups{i},:)];  
      end
    end
    X = X - mean(X);
    beta = Y\X;
    beta = beta./sum(beta);

    BA_unsorted_weighted = sum(BA_unsorted*beta,2);

    % scale estimated BA values by ratio between SD of original and estimated BA values to get the same range
    BA_unsorted_weighted = BA_unsorted_weighted*mean(std(BA_unsorted))./mean(std(BA_unsorted_weighted));
    
    if 0 & D.verbose
      for i=1:numel(beta)
        fprintf('%s:\tweighting %g\n',Names_all{i},beta(i));
      end
    end
  end
    
  % apply trend correction
  BA_unsorted_weighted = apply_trend_correction(BA_unsorted_weighted,age,D);
  BA_weighted = BA_unsorted_weighted(ind_groups);
  MAE_ctl_weighted = mean(abs(BA_unsorted_weighted(D.ind_adjust)));

  fprintf('-----------------------------------------------------------------\n'); 
  fprintf('weighted MAE (testdata): %g\n',MAE_ctl_weighted); 
  fprintf('-----------------------------------------------------------------\n'); 
  
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
      for k = 1:n_groups
        fprintf('%20s\t',deblank(D.name_groups(k,:)));
      end
      fprintf('\n');
    
      for o = 1:n_groups
        fprintf('%20s\t',deblank(D.name_groups(o,:)));
        for p = 1:n_groups
          if P(o,p) <= 0.001 | P(o,p) >= 0.999
            fprintf('%20.7f***\t',P(o,p));
          elseif P(o,p) <= 0.01 | P(o,p) >= 0.99
            fprintf('%20.7f **\t',P(o,p));
          elseif P(o,p) <= 0.05 | P(o,p) >= 0.95
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
end

% show plot for multiple values if defined
if multiple_BA

  % build testmatrix according to number of groups
  switch n_groups
  case {1,2}
    n_test = 1; test_matrix = [1 2];
  case 3
    n_test = 3; test_matrix = [1 2; 1 3; 2 3];
  case 4
    n_test = 6; test_matrix = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
  case 5
    n_test = 10; test_matrix = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5];
  otherwise
    fprintf('Only comparison of max. 5 groups allowed.');
    return
  end
  P_array = zeros(n_age_range,numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array)*numel(threshold_std));
  MAE_array = zeros(n_age_range,numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array)*numel(threshold_std));
  
  name_array = cell(numel(D.smooth_array)*numel(D.seg_array)*numel(D.res_array)*numel(threshold_std),1);
  age_array = cell(n_age_range,1);
  for m=1:n_age_range            
    age_array{m} = sprintf('%g-%g',age_range{m}(1),age_range{m}(2));
  end
  
  count = 0; count0 = 0;
  for i = 1:numel(D.res_array)
    for j = 1:numel(D.smooth_array)
      for o = 1:numel(D.seg_array)
        count = count + 1;  
        for l = 1:numel(threshold_std)
          for o = 1:numel(trend_degree)
        
            count0 = count0 + 1;
            name_array{count0} = [seg '-' D.smooth '-' D.res 'mm-std' num2str(threshold_std(l)) '-trend' num2str(trend_degree(o))];
            for m = 1:n_age_range          
              if exist('P')
                P_array(m,count0) = P_all{count,m,l};
              end
              MAE_array(m,count0) = MAE_ctl{count0};
            end
          end
        end
      end
    end
  end
  
  ind_groups = [];
  for k = 1:n_groups
    ind_groups = [ind_groups; D.ind_groups{k}];
    data_cell{k} = BA_unsorted_weighted(D.ind_groups{k});
    avg_BrainAGE(k) = mean(BA_unsorted_weighted(D.ind_groups{k}));
    age_test(k,1) = mean(D.age_test(D.ind_groups{k}));
    age_test(k,2) = std(D.age_test(D.ind_groups{k}));
    median_BrainAGE(k) = median(BA_unsorted_weighted(D.ind_groups{k}));
  end
  
  if exist('P') && D.verbose
    figure(22)
          
    if n_age_range == 1
      col = num2cell(jet(count0),2);
      plot(P_array(1),'*','Color',col{1});
      hold on
      for i=2:count0
        plot(P_array(i),'*','Color',col{i});
      end
      hold off
    else
      h = plot(P_array);
      set(h, {'color'}, num2cell(jet(count0),2));
    end
      
    legend(name_array)
    set(gca,'XTick',1:n_age_range,'XTickLabel',age_array)
    title('Anova');
  
  end

  figure(23)
        
  if n_age_range == 1
    col = num2cell(jet(count0),2);
    plot(MAE_array(1),'*','Color',col{1});
    hold on
    for i=2:count0
      plot(MAE_array(i),'*','Color',col{i});
    end
    hold off
  else
    h = plot(MAE_array);
    set(h, {'color'}, num2cell(jet(count0),2));
  end
    
  legend(name_array)
  set(gca,'XTick',1:n_age_range,'XTickLabel',age_array)
  title('MAE [years]');
  
  % print BrainAGE of groups
  if D.verbose
    fprintf('%20s\t','Group');
    for k = 1:n_groups
      fprintf('%20s\t',deblank(D.name_groups(k,:)));
    end
  
    fprintf('\n'); fprintf('%20s\t','Mean');   for k = 1:n_groups, fprintf('%20.3f\t',avg_BrainAGE(k)); end
    fprintf('\n'); fprintf('%20s\t','Median'); for k = 1:n_groups, fprintf('%20.3f\t',median_BrainAGE(k)); end
    fprintf('\n');

    figure(24)
    tmp = cat_plot_boxplot(data_cell,opt);

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
function weighted_EstimatedAge = get_EstimatedAge(BA,age)
%-------------------------------------------------------------------------------
n = numel(BA{1});
m = numel(BA);
  
EstimatedAge = zeros(n,m);
for i=1:m
  EstimatedAge(:,i) = BA{i} + age;
end

Beta = EstimatedAge\age;
weighted_EstimatedAge = EstimatedAge*Beta;

%-------------------------------------------------------------------------------
function [BrainAGE, EstimatedAge] = apply_trend_correction(BrainAGE,age,D)
%-------------------------------------------------------------------------------

EstimatedAge = BrainAGE + age;
% apply trend correction for each site separately
for i = 1:max(D.site_adjust)
  ind_site = find(D.site_adjust == i);
  
  % build index for each site
  ind_site_adjust = [];
  for j=1:length(D.ind_adjust)
    if D.site_adjust(D.ind_adjust(j)) == i
      ind_site_adjust = [ind_site_adjust D.ind_adjust(j)];
    end
  end

  % don't correct trend if trend_degree is < 0
  if ~isempty(ind_site_adjust) & D.trend_degree >= 0
  
    if length(ind_site_adjust) < 20
      fprintf('Sample #%d for site-specific trend-correction is too small (n=%d). Use only offset-correction by Median instead.\n',i, length(ind_site_adjust));
      offset = median(BrainAGE(ind_site_adjust));
      BrainAGE(ind_site) = BrainAGE(ind_site) - offset;
      EstimatedAge(ind_site) = EstimatedAge(ind_site) - offset;
    else
      if D.trend_degree > 0
        G = [ones(length(age(ind_site)),1) cg_polynomial(age(ind_site),D.trend_degree)];
        G_indexed = [ones(length(age(ind_site_adjust)),1) cg_polynomial(age(ind_site_adjust),D.trend_degree)];
      else
        G = [ones(length(age(ind_site)),1)];
        G_indexed = [ones(length(age(ind_site_adjust)),1)];
      end
      if D.verbose, fprintf('Remove trend degree %d using %d subjects of site %d.\n',D.trend_degree,length(ind_site_adjust),i); end
      
      % estimate beta only for indexed data (e.g. control subjects)
      beta = pinv(G_indexed)*BrainAGE(ind_site_adjust);
      
      % and remove effects for all data
      Gbeta = G*beta;
      BrainAGE(ind_site) = BrainAGE(ind_site) - Gbeta;
      EstimatedAge(ind_site)  = EstimatedAge(ind_site)  - Gbeta;
    end
  else
    if D.trend_degree >= 0
      fprintf('No subjects found in sample #%d for site-specific trend-correction\n',i);
    end
  end
end

if D.trend_degree >= 0
  avg_BrainAGE = mean(BrainAGE(D.ind_adjust));
  BrainAGE = BrainAGE - avg_BrainAGE;
  EstimatedAge = EstimatedAge - avg_BrainAGE;
end
