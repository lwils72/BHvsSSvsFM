function [Y]=quantile(X,p)
% [Y]=quantile(X,p)
%    - compute the quantiles of a given data set
%    - this is a wrapper for native histogram function
%      using 'cdf' normalization
%    - bootleg work around for not having the stats toolbox
%
%  inputs: X = vector of interest
%          p = vector of quantiles of interest [0,1]
%  outputs: Y is a vector same length as p with the quantiles of X
%

% version history:
% Original: Karen Luttrell, 8 June 2015
%   - bare bones histogram wrapper, column vectors only
%   - automatic 1000 bin histogram calculation
%   - in testing, this gives the correct median to 2 significant figures
%		- good enough for now

  N=1000;
	b=linspace(min(X),max(X),N);
	h=histcounts(X,b,'normalization','cdf');

	for i=1:numel(p)
		[junk,j]=min(abs(h-p(i)));
		Y(i)=b(j);
	end

