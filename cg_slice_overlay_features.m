function cg_slice_overlay_features

addpath ~/spm/slice_overlay

OV.reference_image = '../../avg_wmIXI547.nii'; % anatomical image to underlay
OV.reference_range = [200 1300];                         % intensity range for reference image
OV.opacity = 0.8;                                      % transparence value for overlay (<1)
OV.cmap    = [1-(hot); hot];                                      % colormap for overlay

% name of files
OV.name=str2mat(   'features8mm_pca_p5%',...
                'features3mm_p0.05',...
                'features3mm_p0.1');
                
% range for each file
% use range 0..0 if you want to autoscale range
% If you are using log. scaling, check the highest p-value in the table
% and approximate the range; e.g. for a max. p-value of p=1e-7 and a
% threshold of p<0.001 use a range of [3 7]. Check cg_t2x for details.
% if you are unsure, simply use the autorange option by with a range of [0 0]
% to approximate the range.
% The log-scaled values are calculated by -log10(1-p):
% p-value       -log10(1-P)
%  0.1           1
%  0.05          1.3
%  0.01          2
%  0.001         3
%  0.0001        4

% number of fields in range should be the same as number of files (see above)
% if lower and upper range are equal, then the range will be automatically estimated
OV.range   =[[-1.5 1.5];[-.12 .12];[-.12 .12]];
OV.cmap    = [1-(hot); hot];                                      % colormap for overlay

% log. scaling used (get info from filename)
OV.logP    = zeros(size(OV.name,1));
for i=1:size(OV.name,1)
    if findstr(OV.name(i,:),'logP')
        OV.logP(i) = 1;
    end
end

% selection of slices and orientations
OV.slices_str = str2mat('-50:2:72','-70:2:70','-96:2:62');
OV.transform = str2mat('axial','sagittal','coronal');

cg_slice_overlay(OV)
