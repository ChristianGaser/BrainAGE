function [p,anovatab,stats] = cg_anova1(x,group,displayopt,extra)
%ANOVA1 One-way analysis of variance (ANOVA).
%   ANOVA1 performs a one-way ANOVA for comparing the means of two or more 
%   groups of data. It returns the p-value for the null hypothesis that the
%   means of the groups are equal.
%
%   P = ANOVA1(X,GROUP,DISPLAYOPT)
%   If X is a matrix, ANOVA1 treats each column as a separate group, and
%     determines whether the population means of the columns are equal.
%     This form of ANOVA1 is appropriate when each group has the same
%     number of elements (balanced ANOVA).  GROUP can be a character
%     array or a cell array of strings, with one row per column of
%     X, containing the group names.  Enter an empty array ([]) or
%     omit this argument if you do not want to specify group names.
%   If X is a vector, GROUP must be a categorical variable, vector,
%     string array, or cell array of strings with one group name for
%     each element of X.  X values corresponding to the same value of
%     GROUP are placed in the same group.
%
%   DISPLAYOPT can be 'on' (the default) to display figures
%   containing a standard one-way anova table and a boxplot, or
%   'off' to omit these displays.  Note that the notches in the
%   boxplot provide a test of group medians (see HELP BOXPLOT),
%   and this is not the same as the F test for different means
%   in the anova table.
%
%   [P,ANOVATAB] = ANOVA1(...) returns the ANOVA table values as the
%   cell array ANOVATAB.
%
%   [P,ANOVATAB,STATS] = ANOVA1(...) returns an additional structure
%   of statistics useful for performing a multiple comparison of means
%   with the MULTCOMPARE function.
%
%   See also ANOVA2, ANOVAN, BOXPLOT, MANOVA1, MULTCOMPARE.

%   Reference: Robert V. Hogg, and Johannes Ledolter, Engineering Statistics
%   Macmillan 1987 pp. 205-206.

%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/03/22 04:41:18 $

classical = 1;
nargs = nargin;
if (nargin>0 && strcmp(x,'kruskalwallis'))
   % Called via kruskalwallis function, adjust inputs
   classical = 0;
   if (nargin >= 2), x = group; group = []; end
   if (nargin >= 3), group = displayopt; displayopt = []; end
   if (nargin >= 4), displayopt = extra; end
   nargs = nargs-1;
end

error(nargchk(1,3,nargs,'struct'));

if (nargs < 2), group = []; end
if (nargs < 3), displayopt = 'on'; end
% Note: for backwards compatibility, accept 'nodisplay' for 'off'
willdisplay = ~(strcmp(displayopt,'nodisplay') | strcmp(displayopt,'n') ...
                | strcmp(displayopt,'off'));

% Convert group to cell array from character array, make it a column
if (ischar(group) && ~isempty(group)), group = cellstr(group); end
if (size(group, 1) == 1), group = group'; end

% If X is a matrix with NaNs, convert to vector form.
if (length(x) < numel(x))
   if (any(isnan(x(:))))
      [n,m] = size(x);
      x = x(:);
      gi = reshape(repmat((1:m), n, 1), n*m, 1);
      if isempty(group)     % no group names
         group = gi;
      elseif (size(group,1) == m)
         group = group(gi,:);
      else
         error('stats:anova1:InputSizeMismatch',...
               'X and GROUP must have the same length.');
      end
   end
end

% If X is a matrix and GROUP is strings, use GROUPs as names
if (iscell(group) && (length(x) < numel(x)) ...
                  && (size(x,2) == size(group,1)))
   named = 1;
   gnames = group;
   grouped = 0;
else
   named = 0;
   gnames = [];
   grouped = ~isempty(group);
end

if (grouped)
   % Single data vector and a separate grouping variable
   x = x(:);
   lx = length(x);
   if (lx ~= numel(x))
      error('stats:anova1:VectorRequired','First argument has to be a vector.')
   end
   nonan = ~isnan(x);
   x = x(nonan);

   % Convert group to indices 1,...,g and separate names
   group = group(nonan,:);
   [groupnum, gnames] = grp2idx(group);
   named = 1;

   % Remove NaN values
   nonan = ~isnan(groupnum);
   if (~all(nonan))
      groupnum = groupnum(nonan);
      x = x(nonan);
   end

   lx = length(x);
   xorig = x;                    % use uncentered version to make M
   groupnum = groupnum(:);
   maxi = size(gnames, 1);
   if isa(x,'single')
      xm = zeros(1,maxi,'single');
   else
      xm = zeros(1,maxi);
   end
   countx = xm;
   if (classical)
      mu = mean(x);
      x = x - mu;                % center to improve accuracy
      xr = x;
   else
      [xr,tieadj] = tiedrank(x);
   end
   
   if (willdisplay)
      % Fill M with NaN in case the group sizes vary
      Mrows = 0;
      for j=1:maxi
         Mrows = max(Mrows,sum(groupnum==j));
      end
      M = NaN(Mrows,maxi);
   end
   for j = 1:maxi
      % Get group sizes and means
      k = find(groupnum == j);
      lk = length(k);
      countx(j) = lk;
      xm(j) = mean(xr(k));       % column means

      if (willdisplay)           % fill matrix for boxplot
         M(1:lk,j) = xorig(k);
      end
   end

   gm = mean(xr);                      % grand mean
   df1 = sum(countx>0) - 1;            % Column degrees of freedom
   df2 = lx - df1 - 1;                 % Error degrees of freedom
   xc = xm - gm;                       % centered
   xc(countx==0) = 0;
   RSS = dot(countx, xc.^2);           % Regression Sum of Squares
else
   % Data in matrix form, no separate grouping variable
   [r,c] = size(x);
   lx = r * c;
   if (classical)
      xr = x;
      mu = mean(xr(:));
      xr = xr - mu;           % center to improve accuracy
   else
      [xr,tieadj] = tiedrank(x(:));
      xr = reshape(xr, size(x));
   end
   countx = repmat(r, 1, c);
   xorig = x;                 % save uncentered version for boxplot
   xm = mean(xr);             % column means
   gm = mean(xm);             % grand mean
   df1 = c-1;                 % Column degrees of freedom
   df2 = c*(r-1);             % Error degrees of freedom
   RSS = r*(xm - gm)*(xm-gm)';        % Regression Sum of Squares
end

TSS = (xr(:) - gm)'*(xr(:) - gm);  % Total Sum of Squares
SSE = TSS - RSS;                   % Error Sum of Squares

if (df2 > 0)
   mse = SSE/df2;
else
   mse = NaN;
end

if (classical)
   if (SSE~=0)
      F = (RSS/df1) / mse;
      p = fpval(F,df1,df2);        % Probability of F given equal means.
   elseif (RSS==0)                 % Constant Matrix case.
      F = NaN;
      p = NaN;
   else                            % Perfect fit case.
      F = Inf;
      p = 0;
   end
else
   F = (12 * RSS) / (lx * (lx+1));
   if (tieadj > 0)
      F = F / (1 - 2 * tieadj/(lx^3-lx));
   end
   p = chi2pval(F,df1);
end


Table=zeros(3,5);               %Formatting for ANOVA Table printout
Table(:,1)=[ RSS SSE TSS]';
Table(:,2)=[df1 df2 df1+df2]';
Table(:,3)=[ RSS/df1 mse Inf ]';
Table(:,4)=[ F Inf Inf ]';
Table(:,5)=[ p Inf Inf ]';

colheads = ['Source       ';'         SS  ';'          df ';...
            '       MS    ';'          F  ';'     Prob>F  '];
if (~classical)
   colheads(5,:) = '     Chi-sq  ';
   colheads(6,:) = '  Prob>Chi-sq';
end
rowheads = ['Columns    ';'Error      ';'Total      '];
if (grouped)
   rowheads(1,:) = 'Groups     ';
end

% Create cell array version of table
atab = num2cell(Table);
for i=1:size(atab,1)
   for j=1:size(atab,2)
      if (isinf(atab{i,j}))
         atab{i,j} = [];
      end
   end
end
atab = [cellstr(strjust(rowheads, 'left')), atab];
atab = [cellstr(strjust(colheads, 'left'))'; atab];
if (nargout > 1)
   anovatab = atab;
end

% Create output stats structure if requested, used by MULTCOMPARE
if (nargout > 2)
   if ~isempty(gnames)
      stats.gnames = gnames;
   else
      stats.gnames = strjust(num2str((1:length(xm))'),'left');
   end
   stats.n = countx;
   if (classical)
      stats.source = 'anova1';
      stats.means = xm + mu;
      stats.df = df2;
      stats.s = sqrt(mse);
   else
      stats.source = 'kruskalwallis';
      stats.meanranks = xm;
      stats.sumt = 2 * tieadj;
   end
end

if (~willdisplay), return; end

digits = [-1 -1 0 -1 2 4];
if (classical)
   wtitle = 'One-way ANOVA';
   ttitle = 'ANOVA Table';
else
   wtitle = 'Kruskal-Wallis One-way ANOVA';
   ttitle = 'Kruskal-Wallis ANOVA Table';
end
tblfig = statdisptable(atab, wtitle, ttitle, '', digits);
set(tblfig,'tag','table');

f1 = figure('pos',get(gcf,'pos') + [0,-200,0,0],'tag','boxplot');
ax = axes('Parent',f1);
if (~grouped)
   boxplot(ax,xorig,'notch','on');
else
   boxplot(ax,M,'notch','on');
   h = get(ax,'XLabel');
   set(h,'String','Group Number');
end

% If there are group names, use them
if ~isempty(gnames)
   h = get(ax,'XLabel');
   if (named)
      set(h,'String','');
   end
   set(ax, 'xtick', (1:size(gnames,1)), 'xticklabel', gnames);
end

function [gidx,gnames,glevels] = grp2idx(s)
% GRP2IDX  Create index vector from a grouping variable.
%   [G,GN] = GRP2IDX(S) creates an index vector G from the grouping variable
%   S. S can be a categorical, numeric, or logical vector; a cell vector of
%   strings; or a character matrix with each row representing a group label.
%   The result G is a vector taking integer values from 1 up to the number K
%   of distinct groups. GN is a cell array of strings representing group
%   labels. GN(G) reproduces S (aside from any differences in type).
%
%   Type "help groupingvariable" for more information about grouping
%   variables.
%
%   [G,GN,GL] = GRP2IDX(S) returns a column vector GL representing the
%   group levels. The set of groups and their order in GL and GN are the
%   same, except that GL has the same type as S. If S is a character
%   matrix, GL(G,:) reproduces S, otherwise GL(G) reproduces S.
%
%   GRP2IDX treats NaNs (numeric or logical), empty strings (char or cell
%   array of strings), or <undefined> values (categorical) in S as missing
%   values and returns NaNs in the corresponding rows of G. GN and GL don't
%   include entries for missing values.
%
%   See also GROUPINGVARIABLE, GRPSTATS, GSCATTER.

%   Copyright 1999-2008 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:14:32 $

charFlag = ischar(s);
if charFlag
    charWidth = size(s,2);
    if isempty(s)
        s = cell(0,1);
    else
        s = cellstr(s);
    end
end

if ~isvector(s)
    error('stats:grp2idx:BadGroup',...
          'Grouping variable must be a vector or a character array.');
end
s = s(:);

if isnumeric(s) || islogical(s)
    [glevels,dum,gidx] = unique(s,'first');
    
    % Handle NaN missing values: return NaN group indices
    if ~isempty(glevels) && isnan(glevels(end)) % NaNs are sorted to end
        glevels = glevels(~isnan(glevels));
        gidx(gidx > length(glevels)) = NaN;
    end
    if nargout > 1
        gnames = sprintfc('%g',full(glevels));
    end
    
elseif isa(s,'categorical')
    gidx = double(s); % converts <undefined> to NaN
    if nargout > 1
        gnames = getlabels(s)';
        glevels = getlevels(s)';
    end
    
elseif iscell(s)
    try
        [glevels,ord,gidx] = unique(s,'first');
    catch ME
        if isequal(ME.identifier,'MATLAB:CELL:UNIQUE:InputClass')
            error('stats:grp2idx:GroupTypeIncorrect',...
                  ['A grouping variable must be a categorical, numeric, or logical '....
                   'vector, a cell vector of strings, or a 2D character array.']);
        else
            rethrow(ME);
        end
    end
    
    % Get the "first seen" order of the levels
    [dum,reord] = sort(ord);
    ireord(reord) = 1:length(reord); ireord = ireord(:);
    
    % Handle empty string missing values: return NaN group indices
    if ~isempty(glevels) && strcmp('',glevels(1)) % '' is sorted to beginning
        reord(reord==1) = [];
        ireord = ireord - (ireord > ireord(1));
        ireord(1) = NaN;
    end
    
    % Put the levels back into "first seen" order
    gidx = ireord(gidx(:)); % force a col, even for 0x0
    if nargout > 1
        glevels = glevels(reord(:)); % force a col, even for 0x0
        gnames = glevels;
        if charFlag
            if isempty(s)
                glevels = char(zeros(0,charWidth));
            else
                glevels = char(glevels);
            end
        end
    end
    
else
    error('stats:grp2idx:GroupTypeIncorrect',...
          ['A grouping variable must be a categorical, numeric, or logical '....
           'vector, a cell vector of strings, or a 2D character array.']);
end

function p = fpval(x,df1,df2)
%FPVAL F distribution p-value function.
%   P = FPVAL(X,V1,V2) returns the upper tail of the F cumulative distribution
%   function with V1 and V2 degrees of freedom at the values in X.  If X is
%   the observed value of an F test statistic, then P is its p-value.
%
%   The size of P is the common size of the input arguments.  A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   See also FCDF, FINV.

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.6.

%   Copyright 2009 The MathWorks, Inc. 
%   $Revision: 1.1.8.1 $  $Date: 2010/03/16 00:29:35 $

if nargin < 3, 
    error('stats:fpval:TooFewInputs','Requires three input arguments.'); 
end

p = spm_Fcdf(1./x,df2,df1);