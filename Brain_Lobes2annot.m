function Brain_Lobes2annot
% temporary function to create annot-file from parcellation

matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.data_vol = {'Brain_Lobes.nii,1'};
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.merge_hemi = 1;
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.mesh32k = 1;
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.sample = {'maxabs'};
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.interp = {'linear'};
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.datafieldname = 'maxabs';
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.mapping.rel_equivol_mapping.class = 'GM';
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.mapping.rel_equivol_mapping.startpoint = 0.5;
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.mapping.rel_equivol_mapping.steps = 1;
matlabbatch{1}.spm.tools.cat.stools.vol2surftemp.mapping.rel_equivol_mapping.endpoint = 0.5;
spm_jobman('run',matlabbatch);

catdir = fileparts(which('cat12'));
DK40_lh = fullfile(catdir,'atlases_surfaces_32k','lh.aparc_DK40.freesurfer.annot');
[vert_lh, atlas_lh, colortable_lh] = cat_io_FreeSurfer('read_annotation',DK40_lh);
DK40_rh = fullfile(catdir,'atlases_surfaces_32k','rh.aparc_DK40.freesurfer.annot');
[vert_rh, atlas_rh, colortable_rh] = cat_io_FreeSurfer('read_annotation',DK40_rh);

g = gifti('mesh.maxabs_Brain_Lobes.Template_T1.gii');

label_lh = round(g.cdata(1:32492));
label_lh(atlas_lh == 0) = 15;
label_lh(label_lh < 11) = 0;
regions = unique(label_lh(label_lh > 0));
n_regions = numel(regions);
struct_names = {'Frontal L','Parietal L','Occipital L','Temporal L','Subcortical L'};
table = [[11:15]' zeros(n_regions,3) [11:15]'];
colortable = struct('numEntries',n_regions,'orig_tab','write_annotation_tabl',...
    'struct_names',{struct_names},'table',table);
cat_io_FreeSurfer('write_annotation', 'lh.Brain_Lobes.annot', vert_lh, label_lh, colortable);

label_rh = round(g.cdata(32493:64984));
label_rh(atlas_rh == 0) = 5;
label_rh(label_rh > 5) = 0;
regions = unique(label_rh(label_rh > 0));
n_regions = numel(regions);
struct_names = {'Frontal R','Parietal R','Occipital R','Temporal R','Subcortical R'};
table = [[1:5]' zeros(n_regions,3) [1:5]'];
colortable = struct('numEntries',n_regions,'orig_tab','write_annotation_tabl',...
    'struct_names',{struct_names},'table',table);
cat_io_FreeSurfer('write_annotation', 'rh.Brain_Lobes.annot', vert_rh, label_rh, colortable);

delete('mesh.maxabs_Brain_Lobes.Template_T1.gii')

cat_roi_values2surf('lh.Brain_Lobes.annot',MAE_avg,'mesh.test.gii');
%cat_surf_results('clim', [2 3]);
cat_surf_results('print');
!mv _mesh.test.png MAX_avg.png

return
MAE_gpr = [2.90307 2.56793 2.79641 2.4992 2.34359 2.46631 2.55467 2.78266 2.49035 2.30427];
cat_roi_values2surf('lh.Brain_Lobes.annot',MAE_gpr,'mesh.test.gii');
%cat_surf_results('clim', [2 3]);
cat_surf_results('colorbar');
cat_surf_results('colorbar');
cat_surf_results('print');
!mv _mesh.test.png MAX_gpr.png
