function [Y, label, ind, dim, mappedY, mapping, beta] = cg_input_adjusted_data(V, beta)
% FORMAT [Y, label, ind, dim] = cg_input_adjusted_data(V)

V0.mat = [1 0 0 -90; 0 1 0 -126; 0 0 1 -72; 0 0 0 1];
V0.dim = [181 217 181];

spm2 = 0;
if strcmp(spm('ver'),'SPM2'), spm1 = 1; end

if ~spm2
	brainmask = fullfile(spm('Dir'),'apriori','brainmask.nii');
%	brainmask='F_gt_100.nii'
else
	brainmask = fullfile(spm('Dir'),'apriori','brainmask.mnc');
end
Vm  = spm_vol(brainmask);

% select images for each variate
if nargin == 0
	don = 0;
	n = Inf;
	for i = 1:100,
		if ~spm2
			P = spm_select([0 n],'image',['Select images for sample ' num2str(i) ' (or press done)']);
		else	
			P = spm_get([0 n],'IMAGE',['Select images for sample ' num2str(i) ' (or press done)']);
		end

		if isempty(P), don = 1; break; end;
		V{i} = spm_vol(P);
	end
end

n_samples = length(V);
n_subjects = zeros(n_samples,1);
label = [];
for j=1:n_samples
	n_subjects(j) = size(V{j},1);
	label = [label; j*ones(n_subjects(j),1)];
end
n_subjects_all = sum(n_subjects);

% spatial resolution of images
resolution = spm_input('Spatial resolution','+1','m','1.5mm|2mm|3mm|4mm|8mm|16mm',[1.5 2 3 4 8 16],1);
V0.dim = round(V0.dim/resolution);
V0.mat(1:3,1:3) = resolution*V0.mat(1:3,1:3);
V0.mat(1:3,4) = V0.mat(1:3,4) - [resolution resolution resolution]';

confound = spm_input('Nuisance variables',1,'no|yes',[0 1],1);
if confound
	G = spm_input('Nuisance variables','+1','r',[],[n_subjects_all Inf]);
end

% global scaling
scaling = spm_input('Global scaling','+1','m','None|User specified globals|Compute as mean voxel value',[1 2 3],1);
switch scaling
case 2
	%-User specified globals
	g = spm_input('User specified globals','+1','r',[],[n_subjects_all 1]);
	count = 1;
	for j=1:n_samples     
		for i = 1:n_subjects(j)
			V{j}(i).pinfo(1:2,:) = V{j}(i).pinfo(1:2,:)*100/g(count);
			count = count + 1;
		end
	end
case 3
	%-Compute as mean voxel value (within per image fullmean/8 mask)
	fprintf('Calculating globals\n')
	for j=1:n_samples     
		for i = 1:n_subjects(j)
			g = spm_global(V{j}(i));
			V{j}(i).pinfo(1:2,:) = V{j}(i).pinfo(1:2,:)*100/g;
		end
	end
end

if nargout > 2
	mask = zeros(V0.dim(1:3));
	for sl=1:V0.dim(3)
		% read mask
		M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
		Mm  = V0.mat\Vm.mat\M;
		mask(:,:,sl) = spm_slice_vol(Vm,Mm,V0.dim(1:2),1);
	end
	ind = find(mask > 0.5);
	dim = V0.dim(1:3);
end

Y = [];
Ymean = [];
C = zeros(n_subjects_all);

spm_progress_bar('Init',V0.dim(3),'reading...','planes completed');
for sl=1:V0.dim(3)
	% read mask
	M = spm_matrix([0 0 sl 0 0 0 1 1 1]);
	Mm  = V0.mat\Vm.mat\M;
	M1  = V0.mat\V{1}(1).mat\M;
	mask_slice = spm_slice_vol(Vm,Mm,V0.dim(1:2),1);
	ind0 = find(mask_slice > 0.5);
	clear mask_slice

	% read data inside mask
	if ~isempty(ind0)
		yslice = [];
		for j=1:n_samples
			y = zeros(n_subjects(j), length(ind0));
	    	for i = 1:n_subjects(j)
	   	    	d = spm_slice_vol(V{j}(i),M1,V0.dim(1:2),1);
	       		y(i,:) = d(ind0);
	   		end
			yslice = [yslice; y];
		end
		if confound
			if nargin < 2
				beta_slice = pinv(G)*yslice;
			else
			
				beta_slice = beta{sl};
			end
			if nargout == 7
				beta{sl} = beta_slice;
			end
			for k=1:size(G,2)
				yslice = yslice - G(:,k)*beta_slice(k,:);
			end
		end
		if nargout > 4
			yslicemean = mean(yslice);
	  	for i=1:size(yslice,1)
	  		yslice(i,:) = yslice(i,:) - yslicemean;
			end
			Ymean = [Ymean yslicemean];
			C = C + (1 / size(yslice, 1)) * (yslice * yslice'); 
		end
		Y = [Y yslice];
	end
	spm_progress_bar('Set',sl)
end

spm_progress_bar('Clear')
if nargout > 4
	C = C/V0.dim(3);
	
	% Perform eigendecomposition of C
	C(isnan(C)) = 0;
	C(isinf(C)) = 0;
	
	[M, lambda] = eig(C);

  n_dims = n_subjects_all - 1;    
	% Sort eigenvectors in descending order
	[lambda, ind2] = sort(diag(lambda), 'descend');
	M = M(:,ind2(1:n_dims));
	lambda = lambda(1:n_dims);
	
	% Apply mapping on the data
	if (size(Y, 2) > size(Y, 1))
		M2 = (Y'*M);
		lambda2 = (1 ./ sqrt(size(Y, 1) .* lambda));
		for i=1:n_dims
			M2(:,i) = M2(:,i)*lambda2(i);
		end
		M = M2;
	end
	mappedY = Y * M;
    
	% Store information for out-of-sample extension
	mapping.M = M;
	mapping.lambda = lambda;
	mapping.mean = Ymean;
	
end