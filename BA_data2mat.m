function BA_data2mat(D, fwhm, res, seg)
% BA_DATA2MAT Save spatially registered brain volumes as Matlab matrices.
%
% This function prepares brain imaging data for machine learning analysis by
% saving spatially registered volumes as Matlab data matrices. It is specifically
% designed for applications such as Gaussian Process Regression where the spatial
% structure of the data is not considered. The function applies a mask to the
% volume data to remove non-brain areas, ensuring that only relevant brain
% information is included in the output and resamples and smoothes the data with
% different sizes (i.e. 4/8mm resampling, 4/8mm smoothing).
%
% Usage:
%   BA_data2mat(D)
%
% Input Arguments:
%   D (struct) - A structure containing mandatory and optional fields that
%   describe the data processing parameters.
%
% Mandatory Fields:
%   D.data    - Cell array of strings; paths to data folders to be concatenated.
%   D.release - String; release or version information of the data.
%   D.age     - Numeric array; subjects' ages.
%   D.male    - Logical array; indicates male (true) or not (false) for subjects.
%   D.name    - String; the base name for the output .mat file.
%
% Optional Fields:
%   D.fwhm      - Cell array of numbers; smoothing sizes in millimeters. Default: {4, 8}.
%   D.res       - Cell array of numbers; resampling resolutions in millimeters. Default: {4, 8}.
%   D.seg       - Cell array of strings; segmentations to include. Default: {'rp1', 'rp2'}.
%   D.subfolder - String; subfolder name where segmentations are saved. Default: [seg release].
%   D.add_str   - String; additional string appended to the default subfolder name. Default: ''.
%   D.age_range - Numeric array; specifies the inclusive age range as [minAge, maxAge]. Default: [0, Inf].
%   D.mask_th   - Numeric value; specifies the threshold for the brainmask. Default: 0.5.
%
% Example:
%   % Define processing parameters
%   D.fwhm = {8};
%   D.res = {8};
%   D.seg = {'rp1', 'rp2', 'rp3'};
%   D.age = load('your_age_data.mat');
%   D.data = {'folder1', 'folder2'};
%   D.name = 'test_data';
%   D.release = '_CAT12.9';
%
%   % Execute function
%   BA_data2mat(D);
%
% Output:
%   The function will save segmented, smoothed, and resampled volumes as .mat
%   files in the specified or default directories. For example, using the above
%   parameters, the output files will be:
%     - s8_8mm_rp1_test_data_CAT12.9.mat
%     - s8_8mm_rp2_test_data_CAT12.9.mat
%     - s8_8mm_rp3_test_data_CAT12.9.mat
%
%   The data from the specified folders (e.g., 'folder1' and 'folder2') are
%   concatenated (first 'folder1', then 'folder2') and saved in a single .mat file
%   format for each segmentation.
%
% Note:
%   The function assumes that all provided data folders contain compatible
%   segmentations and that the age and male arrays correspond to the order of
%   subjects in these folders. It is important to verify data consistency and
%   compatibility before running this function to ensure accurate processing and
%   output.
%
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

% add cat12 path if not already done
if ~exist('cat_io_data2mat')
  addpath(fullfile(spm('dir'),'toolbox','cat12'));
end

% set some defaults
if isfield(D,'fwhm'), fwhm_default = D.fwhm;
else, fwhm_default = {4,8}; end

if isfield(D,'res'), res_default = D.res;
else, res_default  = {4,8}; end

if isfield(D,'seg'), seg_default = D.seg;
else, seg_default  = {'rp1','rp2'}; end

if ~isfield(D,'add_str'), D.add_str = ''; end

% go through the different resampling, smoothing sizes and segmentations
if nargin < 2
  for i = 1:numel(fwhm_default)
    for j = 1:numel(res_default)
      for k = 1:numel(seg_default)
        fprintf('\ns%d%s_%dmm:\n',fwhm_default{i},seg_default{k},res_default{j});
        para = sprintf('(%s,%d,%d,''%s'')','D',fwhm_default{i},res_default{j},seg_default{k});
        eval(['BA_data2mat' para])
      end
    end
  end
  return
end

% some fields are always necessary
if ~isfield(D,'data'), error('data field must be defined'); end
if ~isfield(D,'release'), error('release field must be defined'); end
if ~isfield(D,'name'), error('name field must be defined'); end
if isfield(D,'age'), age = D.age; else, error('age field must be defined');end

% convert data field to cell if necessary
if ~iscell(D.data) && ischar(D.data)
  cellstr(D.data);
end

% default subfolder is string such as 'rp1_r1840' based on segmentation and
% release number
if isfield(D,'subfolder'), subfolder = D.subfolder;
else, subfolder = [seg D.release D.add_str]; end

if isfield(D,'male') male = D.male;
else, male = []; end

% find indices within defined age range
if isfield(D,'age_range')
  ind_age = age >= D.age_range(1) & age <= D.age_range(2);
  age = age(ind_age);
  if ~isempty(male)
    male = male(ind_age);
  end
else
  ind_age = [];
end

files = cell(numel(D.data),1);
n = 0;
for i=1:numel(D.data)
  datafolder = fullfile(D.data{i}, subfolder);
  files{i} = spm_select('FPListRec',datafolder,['^' seg]);
  k = (n + 1):(n + size(files{i},1));
  n = n + size(files{i},1);
  
  if isempty(files{i})
    error('No subjects found in %s.\n',datafolder);
  else
    fprintf('%d subjects found in %s.\n',size(files{i},1),datafolder);
  end
  
  % remove files outsie defined age range
  if ~isempty(ind_age)
    files{i} = files{i}(ind_age(k),:);
    fprintf('%d subjects outside defined age range will be removed from %s.\n',sum(~ind_age(k)),datafolder);
    n = n - sum(~ind_age(k));
  end
  
end

if numel(age) ~= n
  fprintf('Only %d of %d values for age found.\n',numel(age),n);
end

if isfield(D,'mask_th')
  name = ['s' num2str(fwhm) seg '_' num2str(res) 'mm_' D.name '_mask' num2str(D.mask_th) D.release '.mat'];
else
  name = ['s' num2str(fwhm) seg '_' num2str(res) 'mm_' D.name D.release '.mat'];
end

data_struct = struct('data',{files},'resolution',res,'fwhm',fwhm,'mask',...
    cat_get_defaults('extopts.brainmask'),'fname',name);
if isfield(D,'mask_th')
  data_struct.mask_th = D.mask_th;
end

% save mat-files using cat_io_data2mat
if isfield(D,'male')
  cat_io_data2mat(data_struct,struct('age',age,'male',male));
else
  cat_io_data2mat(data_struct,struct('age',age));
end
