function show_features(ind,resolution)

load IXI_s8rp3mm_weights

if nargin < 2
  load ind_3mm
  resolution = 3;
end


V0.mat = [1 0 0 -90; 0 1 0 -126; 0 0 1 -72; 0 0 0 1];
V0.dim = [181 217 181];
V0.dim = round(V0.dim/resolution);
V0.mat(1:3,1:3) = resolution*V0.mat(1:3,1:3);
V0.mat(1:3,4) = V0.mat(1:3,4) - [resolution resolution resolution]';
V0.dt = [4 0];
V0.fname = 'features3mm.img';

vol = zeros(V0.dim);
n = length(weig);

vol(ind) = weig(1:n);

spm_write_vol(V0,vol);