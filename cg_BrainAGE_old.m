function [BrainAGE,est_age] = cg_BrainAGE(Y_test,age_test,training_sample,additional_samples,age_range,smooth,res,seg,ind_for_trend_estimation,nuisance,relnumber)
%
% training_sample: 
%           0 - use 410 of 547 from IXI547
%           1 - use IXI547
%          -1 - use only additional sample
%
% additional_samples: 
%           []  - nothing
%           1   - OASIS316
%           2   - NIH394
%           3   - SHIP1000
%           4   - NIH755
%           [1 2] - OASIS316+NIH394
%
% seg: 
%           'rp1'  - GM
%           'rp2'  - GM
%           str2mat('rp1','rp2')  - GM+WM

if nargin < 11, relnumber = '_r432'; end
if nargin < 10, nuisance = []; end
if nargin < 9, ind_for_trend_estimation = 1:length(age_test); end
if nargin < 8, seg = 'rp1'; end
if nargin < 7, res = '8'; end
if nargin < 6, smooth = 's8'; end
if nargin < 5, age_range = [0 Inf]; end
if nargin < 4, additional_samples = []; end
if nargin < 3, training_sample = 1; end

trend_degree = 2;

% add spider paths
addpath /Users/gaser/Desktop/rvm/spider
spider_path;
addpath /Volumes/UltraMax/BrainAGE_core

if training_sample < 0 && isempty(additional_samples)
  error('You have to define at least an additional training sample if no sample is defined.'); 
end

if training_sample >= 0
  load(['/Volumes/UltraMax/IXI-database/BrainAGE_core/' smooth seg(1,:) '_' res 'mm_IXI547' relnumber '.mat'])

  n = length(age);
  [i,ind] = sort(age);
  ind_test = ind(1:4:n);
  ind(1:4:n) = [];
  ind_train = ind;

  agetrain = age(ind_train);
  agetest  = age(ind_test);
  Ytest = Y(ind_test,:);
  Y2    = Y(ind_train,:);
  
  % use additional segmentation if defined
  for i= 2:size(seg,1)
    load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(i,:) '_' res 'mm_IXI547' relnumber '.mat'])
    Ytest = [Ytest Y(ind_test,:)];
    Y2    = [Y2 Y(ind_train,:)];
  end
  
end

switch training_sample
  case 0
    age_train = agetrain;
    IXI_age_test = agetest;
    Y_train = Y2;
    IXI_Y_test = Ytest;
  case 1
    age_train = [agetrain; agetest];
    IXI_age_test = [];
    Y_train = [Y2; Ytest];
    IXI_Y_test = [];
  case -1
    age_train = [];
    IXI_age_test = [];
    Y_train = [];
    IXI_Y_test = [];
end

for i = 1:length(additional_samples)
  switch additional_samples(i)
  case 1
    load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(1,:) '_' res 'mm_OASIS316' relnumber '.mat'])
  case 2
    load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(1,:) '_' res 'mm_NIH394' relnumber '.mat'])
  case 3
    load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(1,:) '_' res 'mm_Ship1000' relnumber '.mat'])
  case 4
    load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(1,:) '_' res 'mm_NIH876' relnumber '.mat'])
  otherwise
    return
 end
 age_train = [age_train; age];
 Y_train = [Y_train; Y];

  % use additional segmentation if defined
  for j= 2:size(seg,1)
    switch additional_samples(i)
    case 1
      load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(j,:) '_' res 'mm_OASIS316' relnumber '.mat'])
    case 2
      load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(j,:) '_' res 'mm_NIH394' relnumber '.mat'])
    case 3
      load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(j,:) '_' res 'mm_Ship1000' relnumber '.mat'])
    case 4
      load(['/Volumes/UltraMax/IXI-database/classification/' smooth seg(j,:) '_' res 'mm_NIH876' relnumber '.mat'])
    otherwise
      return
    end
    Y_train = [Y_train Y];
  end

end

ind_age = find(age_train >= age_range(1) & age_train <= age_range(2));
age_train = age_train(ind_age);
Y_train = Y_train(ind_age,:);

Y = Y_train;
n_subjects = size(Y,1);

% remove age effects from Y to reliably estimate covariance without its influence
Ymean = repmat(mean(Y), [n_subjects 1]);
G = age_train - mean(age_train);
Y = Y - G*(pinv(G)*Y) + Ymean;

% calculate covariance matrix
YpY = (Y*Y')/n_subjects;

% normalize YpY
d      = sqrt(diag(YpY)); % sqrt first to avoid under/overflow
dd     = d*d';
YpY    = YpY./(dd+eps);
t      = find(abs(YpY) > 1); 
YpY(t) = YpY(t)./abs(YpY(t));
YpY(1:n_subjects+1:end) = sign(diag(YpY));

YpYsum = sum(YpY,1);
[iY, jY] = sort(YpYsum, 2, 'descend');

% extract mean correlation for each data set
mean_cov = zeros(n_subjects,1);
for i=1:n_subjects
  % extract row for each subject
  cov0 = YpY(i,:);

  % remove cov with its own
  cov0(i) = [];
  mean_cov(i) = mean(cov0);
end

% threshold for std
threshold_std = 10;

% sort files
[mean_cov_sorted, ind_sorted] = sort(mean_cov,'descend');
threshold_cov = mean(mean_cov) - threshold_std*std(mean_cov);
n_thresholded = min(find(mean_cov_sorted < threshold_cov));

ind_removed = find(mean_cov < threshold_cov);
fprintf('%d subjects removed because their mean covariance was below %g std.\n',length(ind_removed),threshold_std);

age_train(ind_removed) = [];
Y_train(ind_removed,:) = [];

fprintf('%d subjects used for training (age %3.1f..%3.1f years)\n',length(age_train),min(age_train),max(age_train));


%--- KERNEL -> RVR
s = relvm_r(kernel('poly',1));
% get rid of excessive output
s.algorithm.verbosity = 0;

%----- range 0..1
mn = min(Y_train(:));
mx = max(Y_train(:));
Y_train = (Y_train-mn)/(mx-mn);

mn = min(Y_test(:));
mx = max(Y_test(:));
Y_test = (Y_test-mn)/(mx-mn);

if ~training_sample
 mn = min(IXI_Y_test(:));
 mx = max(IXI_Y_test(:)); 
 Y_test_IXI = (IXI_Y_test-mn)/(mx-mn);
end

%---PCA---
[mapped_train, mapping] = pca(Y_train, size(Y_train, 1)-1);
mapped_test = (Y_test - repmat(mapping.mean, [size(Y_test, 1) 1])) * mapping.M;
if ~training_sample
 mapped_test_IXI = (Y_test_IXI - repmat(mapping.mean, [size(Y_test_IXI, 1) 1])) * mapping.M;
end

% range 0..1
mn = min(mapped_train(:));
mx = max(mapped_train(:));
mapped_train = (mapped_train - mn)/(mx - mn);
mapped_test = (mapped_test - mn)/(mx - mn);
if ~training_sample
 mapped_test_IXI = (mapped_test_IXI - mn)/(mx - mn);
end

d = data(mapped_train, age_train);
t1 = data(mapped_test, age_test);
if ~training_sample
 t_IXI = data(mapped_test_IXI, IXI_age_test);
end

%---RVR---
[est,model] = train(s,d);
% Test
pred1 = test(model,t1);
est_age = pred1.X;
BrainAGE = est_age-age_test;

if ~isempty(nuisance)
  G = [ones(length(age_test),1) nuisance];
  fprintf('Remove effect of nuisance parameter\n');

  % estimate beta
  beta = pinv(G)*BrainAGE;

  % and remove effects
  BrainAGE = BrainAGE - G*beta;
  est_age = est_age - G*beta;

end

G = [ones(length(age_test),1) cg_polynomial(age_test,trend_degree)];
G_indexed = [ones(length(age_test(ind_for_trend_estimation)),1) cg_polynomial(age_test(ind_for_trend_estimation),trend_degree)];
fprintf('Remove trend degree %d\n',trend_degree);

% estimate beta only for indexed data (e.g. control subjects)
beta = pinv(G_indexed)*BrainAGE(ind_for_trend_estimation);

% and remove effects for all data
BrainAGE = BrainAGE - G*beta;

est_age = est_age - G*beta;

if ~training_sample
 pred_IXI = test(model,t_IXI);
 est_age_IXI = pred_IXI.X;
 testage_IXI = pred_IXI.Y;
 diff_IXI = mean(testage_IXI)-mean(est_age_IXI);
 est_age_IXI = est_age_IXI+diff_IXI;
 cc = corrcoef(est_age_IXI,testage_IXI);
 fprintf('\n corr (IXItest): %1.3f, MAE (IXItest): %2.3f\n',cc(1,2),mean(abs(est_age_IXI-testage_IXI)));
end

cc = corrcoef(est_age,age_test);
fprintf('\n corr (testdata): %1.3f, MAE (testdata): %2.3f\n\n',cc(1,2),mean(abs(est_age(ind_for_trend_estimation)-age_test(ind_for_trend_estimation))));

return