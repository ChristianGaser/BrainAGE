function [mappedX, mapping] = cg_pca(X, no_dims, method)
%PCA Perform the PCA algorithm
%
%   [mappedX, mapping] = pca(X, no_dims, method)
%
% The function runs PCA on a set of datapoints X. The variable
% no_dims sets the number of dimensions of the feature points in the 
% embedded feature space (no_dims >= 1, default = 2). 
% For no_dims, you can also specify a number between 0 and 1, determining 
% the amount of variance you want to retain in the PCA step.
% The function returns the locations of the embedded trainingdata in 
% mappedX. Furthermore, it returns information on the mapping in mapping.
%
% Method can be:
%. 'eig' Eigenvalue Decomposition of the covariance matrix (faster but less accurate, for compatibiliy)
%. 'svd' Singular Value Decomposition of X (the default)
%
% This file is part of the Matlab Toolbox for Dimensionality Reduction v0.3b.
% The toolbox can be obtained from http://www.cs.unimaas.nl/l.vandermaaten
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten
% Maastricht University, 2007
%
% Modifed version pca.m from Laurens van der Maaten to also allow svd which is
% more accurate and ensures compatibility between systems
% ______________________________________________________________________
%
% Christian Gaser
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id$

if ~exist('no_dims', 'var')
		no_dims = 2;
end

if ~exist('method', 'var')
		method = 'eig';
end

% Make sure data is zero mean
mapping.mean = mean(X, 1);
X = X - repmat(mapping.mean, [size(X, 1) 1]);

switch method
	case 'svd'

    X(isnan(X)) = 0;
    X(isinf(X)) = 0;

    [U,sigma,M] = svd(X,'econ');

		if no_dims == 1     % sigma might have only 1 row
				sigma = sigma(1);
		else
				sigma = diag(sigma);
    end
    
    % use largest possible # of dimensions
    if no_dims > size(M, 2)
        no_dims = size(M, 2);
        fprintf('Target dimensionality reduced to %d.\n',no_dims);
    end
    
    lambda = sigma;    
    M = M(:,1:no_dims);
    lambda = lambda(1:no_dims);

    mappedX = X * M;
        
	case 'eig'
    % Compute covariance matrix
    if size(X, 2) < size(X, 1)
        C = cov(X);
    else
        C = (1 / size(X, 1)) * (X * X');        % if N>D, we better use this matrix for the eigendecomposition
    end
    
    % Perform eigendecomposition of C
    C(isnan(C)) = 0;
    C(isinf(C)) = 0;
    [M, lambda] = eig(C);
            
    % Sort eigenvectors in descending order
    [lambda, ind] = sort(diag(lambda), 'descend');
    if no_dims > size(M, 2)
        no_dims = size(M, 2);
        fprintf('Target dimensionality reduced to %d.\n',no_dims);
    end
    M = M(:,ind(1:no_dims));
    lambda = lambda(1:no_dims);
    
		% Force small negative eigenvalues to zero because of rounding error
		lambda((lambda<0)&(abs(lambda)<(eps(lambda(1))*length(lambda)))) = 0;
	
		% Check if eigvalues are all positive
		if any(lambda<0)
			error(message('PCA:CovNotPositiveSemiDefinite'));
		end

    % Apply mapping on the data
    if ~(size(X, 2) < size(X, 1))
        M = (X' * M) .* repmat((1 ./ sqrt(size(X, 1) .* lambda))', [size(X, 2) 1]);     % normalize in order to get eigenvectors of covariance matrix
    end
    mappedX = X * M;
    
  end
  
% Store information for out-of-sample extension
mapping.M = M;
mapping.lambda = lambda;
    