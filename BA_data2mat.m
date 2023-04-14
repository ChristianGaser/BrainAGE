function BA_data2mat(D,fwhm,res,seg)
% Save spatially registered volume as Matlab data matrix for further use with 
% machine learning tools such as Gaussian Process Rgression. Spatial structure 
% of the data is not considered. Volume data will be masked to remove 
% non-brain areas.
% 
% FORMAT BA_data2mat(D)
%
% Mandatory arguments:
% D.data      - cell array of data folder(s) that will be concatenated
% D.release   - confounds data to be removed
% D.age       - filename for saving mat-file
% D.name      - output directory of saved mat-file
%
% Optional arguments:
% D.fwhm      - cell array of smoothing sizes (default: {4,8})
% D.res       - cell array of resampling resolutions (default: {4,8})
% D.seg       - cell array of segmentations  (default: {'rp1','rp2'})
% D.subfolder - subfolder where segmentations are saved  (default: [seg release])
%
% Example:
% D.fwhm = {8};
% D.res  = {8};
% D.seg  = {'rp1','rp2','rp3'};
% D.age  = load(your_age_data);
% D.data = {folder1,folder2};
% D.name = 'test_data';
% D.release ='_CAT12.8';
% BA_data2mat(D);
%
% The segmentations rp1-rp3 from folder1 and folder2 will be saved after smoothing
% with 8mm FWHM and resampling to 8mm as mat-files:
% s8_8mm_rp1_test_data_CAT12.8.mat
% s8_8mm_rp2_test_data_CAT12.8.mat
% s8_8mm_rp3_test_data_CAT12.8.mat
% The data from folder1 and folder 2 are concatenated (first folder1, followed 
% by folder2) and saved in single format.
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

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

% go through the different resampling, smoothing sizes and segmentations
if nargin < 2
  for i = 1:numel(fwhm_default)
    for j = 1:numel(res_default)
      for k = 1:numel(seg_default)
        fprintf('\ns%d%s_%dmm: ',fwhm_default{i},seg_default{k},res_default{j});
        para = sprintf('(%s,%d,%d,''%s'')','D',fwhm_default{i},res_default{j},seg_default{k});
        eval(['save_data_mat' para])
      end
    end
  end
  return
end

% some fields are always necessary
if ~isfield(D,'data'), error('data field must be defined'); end
if ~isfield(D,'release'), error('release field must be defined'); end
if ~isfield(D,'name'), error('name field must be defined'); end
if isfield(D,'age'), age = D.age;
else, error('age field must be defined');end

% convert data field to cell if necessary
if ~iscell(D.data) && ischar(D.data)
  cellstr(D.data);
end

% default subfolder is string such as 'rp1_r1840' based on segmentation and
% release number
if isfield(D,'subfolder'), subfolder = D.subfolder;
else, subfolder = [seg D.release]; end

if isfield(D,'male'), male = D.male;
else, male = []; end

files = cell(numel(D.data),1);
n = 0;
for i=1:numel(D.data)
  datafolder = fullfile(D.data{i}, subfolder);
  files{i} = spm_select('FPList',datafolder,['^' seg]);
  n = n + size(files{i},1);
  if isempty(files{i})
    error('No subjects found in %s\n',datafolder);
  else
    fprintf('%d subjects found\n',size(files{i},1));
  end
end

if numel(age) ~= n
  fprintf('Only %d of %d values for age found.\n',numel(age),n);
end

name = ['s' num2str(fwhm) seg '_' num2str(res) 'mm_' D.name D.release];

cat_io_data2mat(struct('data',{files},'resolution',res,'fwhm',fwhm,'mask',cat_get_defaults('extopts.brainmask'),...
   'fname',name),struct('age',age,'male',male));
