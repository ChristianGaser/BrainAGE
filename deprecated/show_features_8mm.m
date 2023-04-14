function show_features

load IXI_pca_vweights
load ind_8mm

resolution = 8;

V0.mat = [1 0 0 -90; 0 1 0 -126; 0 0 1 -72; 0 0 0 1];
V0.dim = [181 217 181];
V0.dim = round(V0.dim/resolution);
V0.mat(1:3,1:3) = resolution*V0.mat(1:3,1:3);
V0.mat(1:3,4) = V0.mat(1:3,4) - [resolution resolution resolution]';
V0.dt = [4 0];
V0.fname = 'features8mm_pca.img';

vol = zeros(V0.dim);
n = length(vweig);

vol(ind) = vweig(1:n);

spm_write_vol(V0,vol);