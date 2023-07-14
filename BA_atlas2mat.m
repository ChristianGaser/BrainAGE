function BA_atlas2mat
% Save spatially registered atlas as Matlab data matrix for further use with 
% BrainAGE tools to estimate regional BrainAGE values. 
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________

% add cat12 path if not already done
if ~exist('cat_defaults')
  addpath(fullfile(spm('dir'),'toolbox','cat12'));
end

% use lobe atlas from Toro
V = spm_vol('Brain_Lobes.nii');
[pth,nam,ext] = spm_fileparts(V.fname);

% assume that we use the standard brainmask from CAT12
brainmask = cat_get_defaults('extopts.brainmask');
Vm  = spm_vol(brainmask{1});

% go through 4 and 8mm spatial resolution
for resolution = [4 8]
  Vres.mat = [1 0 0 -90; 0 1 0 -126; 0 0 1 -72; 0 0 0 1];
  Vres.dim = [181 217 181];
  Vres.dim = round(Vres.dim/resolution);
  Vres.mat(1:3,1:3) = resolution*Vres.mat(1:3,1:3);
  Vres.mat(1:3,4) = Vres.mat(1:3,4) - [resolution resolution resolution]';
  dim = Vres.dim(1:3);
  
  mask_ind = zeros(Vres.dim,'logical');
  
  atlas = [];
  for sl=1:Vres.dim(3)
    % read mask
    M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
    
    Mm  = Vres.mat\Vm.mat\M;
    mask_slice = spm_slice_vol(Vm,Mm,Vres.dim(1:2),1);
    mask_ind(:,:,sl) = mask_slice > 0.5;
    M1 = Vres.mat\V.mat\M;
    
    % use NN interpolation which is not optimal, but prevents wrong interpolation
    % between neighbouring labels
    d = round(spm_slice_vol(V,M1,Vres.dim(1:2),0));
    atlas = [atlas; d(mask_ind(:,:,sl))];
  end
  
  out_name = fullfile(pth,[nam '_' num2str(resolution) 'mm.mat']);
  save(out_name, 'atlas')
end

