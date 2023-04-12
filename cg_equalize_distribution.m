function [sample_sel, c] = cg_equalize_distribution(sample_ref, sample_src, opts)
% [sample_sel, c] = cg_equalize_distribution(sample_ref, sample_src, opts)
% Equals distribution of sample_src w.r.t. sample_ref by creating a cost
% matrix based on the euclidean norm between sample_src and sample_ref.
% This cost matrix is used to solve the assignment problem based on the 
% Hungarian method.
% To speed up the very slow search for huge data sets we internally use rounded
% integer values (after multiplying with 1000) for the cost matrix which is much
% faster and should not influence the results. 
%
% sample_ref  - sample that is used as reference distribution. This can be a
%               matrix of size m x n where m are the dimensions and n the 
%               number of entries.
% sample_src  - sample that is used as input (source) with size m x n1 where
%               m are the dimensions and n1 the number of entries
%
% opts field:
% opts.weight - vector of size m that allows to weight the cost function
% opts.range  - matrix m x 2 which defines the range of sample_src
% opts.tol    - vector of size m that defines tolerance between mean value of
%               sample_ref and the mean value of output sample_sel
%
% Output
% -------
% sample_sel  - selected sample within opts.range and within opts.tol 
%               compared to the mean 
% c           - vector of selected sample values
%_______________________________________________________________________
% Christian Gaser
% $Id: cg_equalize_distribution.m $

% Uses the function hungarian from Niclas Borlin for Hungarian Method
% hungarian.m v1.0  96-06-14

if nargin < 3 || ~isfield(opts,'tol')
  fprintf('At least opts.tol should be defined\n\n');
  help(mfilename)
  if nargout, sample_sel = []; c = []; end
  return
end

[n_ref,m0] = size(sample_ref);
[n_src,m]  = size(sample_src);

if m0 > n_ref
  sample_ref = sample_ref';
  [n_ref,m0] = size(sample_ref);
end

if m > n_src
  trans_src = true;
  sample_src = sample_src';
  [n_src,m]  = size(sample_src);
else
  trans_src = false;
end

if n_ref > n_src, error('Size of reference sample should be smaller than your source sample'); end

if ~isfield(opts,'weight'), opts.weight = ones(m,1); end
if ~isfield(opts,'tol'),    opts.tol = inf(m,1); end
if ~isfield(opts,'range'),  opts.range = repmat([-Inf; Inf],1,m); end

if m ~= size(opts.range,2)
  opts.range = opts.range';
end

if m ~=  m0, error('Size mismatch between reference and source sample.'); end

m = numel(opts.weight);
if m ~=  m0, error('Size mismatch. Weight must be defined with %d entries.',m0'); end
m = numel(opts.tol);
if m ~=  m0, error('Size mismatch. Tolerance must be defined with %d entries.',m0); end

cost = zeros(n_src);
scl = opts.weight./max([sample_ref; sample_src]);

for j = 1:m
  ind = sample_src(:,j) < opts.range(1,j) | sample_src(:,j) > opts.range(2,j);
  sample_src(ind,j) = Inf;
end
  
for i = 1:n_src
  for j = 1:m
    if isfinite(sample_src(i,j))
      cost(i,1:n_ref) = cost(i,1:n_ref) + (scl(j)*(sample_ref(:,j) - sample_src(i,j)).^2)';  
    else
      cost(i,1:n_ref) = 1e9;
    end
  end
  cost(i,1:n_ref) = cost(i,1:n_ref).^0.5;
end

% Use rounded values for cost matrix A and scale with 1000 for rounding issues
% Using integers is much faster....
c = hungarian(round(1000*cost));

n = 1;
while 1
  
  ind = n*n_ref+1:n_src;
  cind = c(ind);

  sample_rem  = sample_src(cind,:);
  n_rem = size(sample_rem,1);
  if n_rem < n_ref, break; end
  
  cost_rem = zeros(n_rem);
  for i = 1:n_rem
    for j = 1:m
      if isfinite(sample_rem(i,j))
       cost_rem(i,1:n_ref) = cost_rem(i,1:n_ref) + (scl(j)*(sample_ref(:,j) - sample_rem(i,j)).^2)';    
      else
        cost_rem(i,1:n_ref) = 1e9;
      end
    end
    cost_rem(i,1:n_ref) = cost_rem(i,1:n_ref).^0.5;
  end

  % Use rounded values for cost matrix A and scale with 1000 for rounding issues
  % Using integers is much faster....
  c_rem = hungarian(round(1000*cost_rem));
  ind_rem = n*n_ref+1:n_rem;
  c_rem_ind = c_rem(ind_rem);
  % sort columns without zero entries to rank the remaining costs in the
  % index
  [~,xi] = sort(min(cost_rem(c_rem_ind,1:n_ref),[],2));
  c_rem(ind_rem) = c_rem_ind(xi);

  c(ind) = cind(c_rem);
  n = n + 1;
end

% go through remaining data
if numel(cind) < n_ref
  % select only those columns that exceed n_ref
  % sort columns without zero entries to rank the remaining costs in the
  % index
  [~,xi] = sort(min(cost(cind,1:n_ref),[],2));
  c(ind) = cind(xi);
end

sample_sel = sample_src(c,:);

% find selected sample inside tolerance
for i = n_ref:n_src
  for j = 1:m
    if abs(mean(sample_sel(1:i,j)) - mean(sample_ref(:,j))) > opts.tol(j) && opts.weight(j)
      c = c(1:i-1);
      sample_sel = sample_sel(1:i-1,:);
      if trans_src, sample_sel = sample_sel'; end
      fprintf('Sample    \tSize\tMean\tSD\n');
      for k = 1:m
        fprintf('Reference \t%d\t%3.2f\t%3.2f\n',size(sample_ref,1),mean(sample_ref(:,k)),std(sample_ref(:,k)));
        fprintf('Source    \t%d\t%3.2f\t%3.2f\t',size(sample_sel,1),mean(sample_sel(:,k)),std(sample_sel(:,k)));
        fprintf('\n');
      end
      fprintf('\n');
      return;
    end
  end
end
if trans_src, sample_sel = sample_sel'; end

function [C,T] = hungarian(A)
%HUNGARIAN Solve the Assignment problem using the Hungarian method.
%
%[C,T] = hungarian(A)
%A - a square cost matrix.
%C - the optimal assignment.
%T - the cost of the optimal assignment.
%s.t. T = trace(A(C,:)) is minimized over all possible assignments.

% Adapted from the FORTRAN IV code in Carpaneto and Toth, "Algorithm 548:
% Solution of the assignment problem [H]", ACM Transactions on
% Mathematical Software, 6(1):104-111, 1980.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.
%                 Department of Computing Science, UmeÃ¥ University,
%                 Sweden. 
%                 All standard disclaimers apply.

% A substantial effort was put into this code. If you use it for a
% publication or otherwise, please include an acknowledgement or at least
% notify me by email. /Niclas

[m,n] = size(A);

if (m~= n)
    error('HUNGARIAN: Cost matrix must be square!');
end

% Save original cost matrix.
if nargout > 1
    orig = A;
end

% Reduce matrix
A = hminired(A);

% Do an initial assignment.
[A,C,U] = hminiass(A);

% Repeat while we have unassigned rows.
while (U(n+1))
    % Start with no path, no unchecked zeros, and no unexplored rows.
    LR = zeros(1,n);
    LC = zeros(1,n);
    CH = zeros(1,n);
    RH = [zeros(1,n) -1];
    
    % No labelled columns.
    SLC = [];
    
    % Start path in first unassigned row.
    r = U(n+1);
    % Mark row with end-of-path label.
    LR(r) = -1;
    % Insert row first in labelled row set.
    SLR = r;
    
    % Repeat until we manage to find an assignable zero.
    while (1)
        % If there are free zeros in row r
        if (A(r,n+1)~= 0)
            % ...get column of first free zero.
            l = -A(r,n+1);
            
            % If there are more free zeros in row r and row r in not
            % yet marked as unexplored..
            if (A(r,l)~= 0 && RH(r) ==0)
                % Insert row r first in unexplored list.
                RH(r) = RH(n+1);
                RH(n+1) = r;
                
                % Mark in which column the next unexplored zero in this row
                % is.
                CH(r) = -A(r,l);
            end
        else
            % If all rows are explored..
            if (RH(n+1) <= 0)
                % Reduce matrix.
                [A,CH,RH] = hmreduce(A,CH,RH,LC,LR,SLC,SLR);
            end
            
            % Re-start with first unexplored row.
            r = RH(n+1);
            % Get column of next free zero in row r.
            l = CH(r);
            % Advance "column of next free zero".
            CH(r) = -A(r,l);
            % If this zero is last in the list..
            if (A(r,l) ==0)
                % ...remove row r from unexplored list.
                RH(n+1) = RH(r);
                RH(r) = 0;
            end
        end
        
        % While the column l is labelled, i.e. in path.
        while (LC(l)~= 0)
            % If row r is explored..
            if (RH(r) ==0)
                % If all rows are explored..
                if (RH(n+1) <= 0)
                    % Reduce cost matrix.
                    [A,CH,RH] = hmreduce(A,CH,RH,LC,LR,SLC,SLR);
                end
                
                % Re-start with first unexplored row.
                r = RH(n+1);
            end
            
            % Get column of next free zero in row r.
            l = CH(r);
            
            % Advance "column of next free zero".
            CH(r) = -A(r,l);
            
            % If this zero is last in list..
            if(A(r,l) ==0)
                % ...remove row r from unexplored list.
                RH(n+1) = RH(r);
                RH(r) = 0;
            end
        end
        
        % If the column found is unassigned..
        if (C(l) ==0)
            % Flip all zeros along the path in LR,LC.
            [A,C,U] = hmflip(A,C,LC,LR,U,l,r);
            % ...and exit to continue with next unassigned row.
            break;
        else
            % ...else add zero to path.
            
            % Label column l with row r.
            LC(l) = r;
            
            % Add l to the set of labelled columns.
            SLC = [SLC l];
            
            % Continue with the row assigned to column l.
            r = C(l);
            
            % Label row r with column l.
            LR(r) = l;
            
            % Add r to the set of labelled rows.
            SLR = [SLR r];
        end
    end
end

% Calculate the total cost.
if nargout > 1
    T = sum(orig(logical(sparse(C,1:size(orig,2),1))));
end

function A = hminired(A)
%HMINIRED Initial reduction of cost matrix for the Hungarian method.
%
%B = assredin(A)
%A - the unreduced cost matris.
%B - the reduced cost matrix with linked zeros in each row.

% v1.0  96-06-13. Niclas Borlin, niclas@cs.umu.se.

[m,n] = size(A);

% Subtract column-minimum values from each column.
colMin = min(A);

A = A-colMin(ones(n,1),:);

% Subtract row-minimum values from each row.
rowMin = min(A')';
A = A-rowMin(:,ones(1,n));

% Get positions of all zeros.
[i,j] = find(A ==0);

% Extend A to give room for row zero list header column.
A(1,n+1) = 0;
for k = 1:n
    % Get all column in this row. 
    cols = j(k ==i)';
    % Insert pointers in matrix.
    A(k,[n+1 cols]) = [-cols 0];
end


function [A,C,U] = hminiass(A)
%HMINIASS Initial assignment of the Hungarian method.
%
%[B,C,U] = hminiass(A)
%A - the reduced cost matrix.
%B - the reduced cost matrix, with assigned zeros removed from lists.
%C - a vector. C(J) = I means row I is assigned to column J,
%              i.e. there is an assigned zero in position I,J.
%U - a vector with a linked list of unassigned rows.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

[n,np1] = size(A);

% Initalize return vectors.
C = zeros(1,n);
U = zeros(1,n+1);

% Initialize last/next zero "pointers".
LZ = zeros(1,n);
NZ = zeros(1,n);

for i = 1:n
    % Set j to first unassigned zero in row i.
  lj = n+1;
  j = -A(i,lj);

    % Repeat until we have no more zeros (j ==0) or we find a zero
  % in an unassigned column (c(j) ==0).
    
  while (C(j)~= 0)
    % Advance lj and j in zero list.
    lj = j;
    j = -A(i,lj);
  
    % Stop if we hit end of list.
    if (j ==0)
      break;
    end
  end

  if (j~= 0)
    % We found a zero in an unassigned column.
    
    % Assign row i to column j.
    C(j) = i;
    
    % Remove A(i,j) from unassigned zero list.
    A(i,lj) = A(i,j);

    % Update next/last unassigned zero pointers.
    NZ(i) = -A(i,j);
    LZ(i) = lj;

    % Indicate A(i,j) is an assigned zero.
    A(i,j) = 0;
  else
    % We found no zero in an unassigned column.

    % Check all zeros in this row.

    lj = n+1;
    j = -A(i,lj);
    
    % Check all zeros in this row for a suitable zero in another row.
    while (j~= 0)
      % Check the in the row assigned to this column.
      r = C(j);
      
      % Pick up last/next pointers.
      lm = LZ(r);
      m = NZ(r);
      
      % Check all unchecked zeros in free list of this row.
      while (m~= 0)
        % Stop if we find an unassigned column.
        if (C(m) ==0)
          break;
        end
        
        % Advance one step in list.
        lm = m;
        m = -A(r,lm);
      end
      
      if (m ==0)
        % We failed on row r. Continue with next zero on row i.
        lj = j;
        j = -A(i,lj);
      else
        % We found a zero in an unassigned column.
      
        % Replace zero at (r,m) in unassigned list with zero at (r,j)
        A(r,lm) = -j;
        A(r,j) = A(r,m);
      
        % Update last/next pointers in row r.
        NZ(r) = -A(r,m);
        LZ(r) = j;
      
        % Mark A(r,m) as an assigned zero in the matrix . . .
        A(r,m) = 0;
      
        % ...and in the assignment vector.
        C(m) = r;
      
        % Remove A(i,j) from unassigned list.
        A(i,lj) = A(i,j);
      
        % Update last/next pointers in row r.
        NZ(i) = -A(i,j);
        LZ(i) = lj;
      
        % Mark A(r,m) as an assigned zero in the matrix . . .
        A(i,j) = 0;
      
        % ...and in the assignment vector.
        C(j) = i;
        
        % Stop search.
        break;
      end
    end
  end
end

% Create vector with list of unassigned rows.

% Mark all rows have assignment.
r = zeros(1,n);
rows = C(C~= 0);
r(rows) = rows;
empty = find(r ==0);

% Create vector with linked list of unassigned rows.
U = zeros(1,n+1);
U([n+1 empty]) = [empty 0];


function [A,C,U] = hmflip(A,C,LC,LR,U,l,r)
%HMFLIP Flip assignment state of all zeros along a path.
%
%[A,C,U] = hmflip(A,C,LC,LR,U,l,r)
%Input:
%A   - the cost matrix.
%C   - the assignment vector.
%LC  - the column label vector.
%LR  - the row label vector.
%U   - the 
%r,l - position of last zero in path.
%Output:
%A   - updated cost matrix.
%C   - updated assignment vector.
%U   - updated unassigned row list vector.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n = size(A,1);

while (1)
    % Move assignment in column l to row r.
    C(l) = r;
    
    % Find zero to be removed from zero list..
    
    % Find zero before this.
    m = find(A(r,:) ==-l);
    
    % Link past this zero.
    A(r,m) = A(r,l);
    
    A(r,l) = 0;
    
    % If this was the first zero of the path..
    if (LR(r)<0)
        ...remove row from unassigned row list and return.
        U(n+1) = U(r);
        U(r) = 0;
        return;
    else
        
        % Move back in this row along the path and get column of next zero.
        l = LR(r);
        
        % Insert zero at (r,l) first in zero list.
        A(r,l) = A(r,n+1);
        A(r,n+1) = -l;
        
        % Continue back along the column to get row of next zero in path.
        r = LC(l);
    end
end

function [A,CH,RH] = hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%HMREDUCE Reduce parts of cost matrix in the Hungerian method.
%
%[A,CH,RH] = hmreduce(A,CH,RH,LC,LR,SLC,SLR)
%Input:
%A   - Cost matrix.
%CH  - vector of column of 'next zeros' in each row.
%RH  - vector with list of unexplored rows.
%LC  - column labels.
%RC  - row labels.
%SLC - set of column labels.
%SLR - set of row labels.
%
%Output:
%A   - Reduced cost matrix.
%CH  - Updated vector of 'next zeros' in each row.
%RH  - Updated vector of unexplored rows.

% v1.0  96-06-14. Niclas Borlin, niclas@cs.umu.se.

n = size(A,1);

% Find which rows are covered, i.e. unlabelled.
coveredRows = LR ==0;

% Find which columns are covered, i.e. labelled.
coveredCols = LC~= 0;

r = find(~coveredRows);
c = find(~coveredCols);

% Get minimum of uncovered elements.
m = min(min(A(r,c)));

% Subtract minimum from all uncovered elements.
A(r,c) = A(r,c)-m;

% Check alluncovered rows in path order..
for i = SLR
    % Check all uncovered columns..
    for j = c
        % If this is a (new) zero..
        if (A(i,j) ==0)
            % If the row is not in unexplored list..
            if (RH(i) ==0)
                % ...insert it first in unexplored list.
                RH(i) = RH(n+1);
                RH(n+1) = i;
                % Mark this zero as "next free" in this row.
                CH(i) = j;
            end
            % Find last unassigned zero on row I.
            row = A(i,:);
            colsInList = -row(row<0);
            if (isempty(colsInList))
                % No zeros in the list.
                l = n+1;
            else
                l = colsInList(row(colsInList) ==0);
            end
            % Append this zero to end of list.
            A(i,l) = -j;
        end
    end
end

% Add minimum to all doubly covered elements.
r = find(coveredRows);
c = find(coveredCols);

% Take care of the zeros we will remove.
[i,j] = find(A(r,c) <= 0);

i = r(i);
j = c(j);

for k = 1:length(i)
    % Find zero before this in this row.
    lj = A(i(k),:) ==-j(k);
    % Link past it.
    A(i(k),lj) = A(i(k),j(k));
    % Mark it as assigned.
    A(i(k),j(k)) = 0;
end

A(r,c) = A(r,c)+m;