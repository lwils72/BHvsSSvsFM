function x = xsym(m)
%XSYM   Symmetric x-axis about zero
%   XSYM(M) changes x axis limits to [-M,M].
%   XSYM, by itself, chooses M to be the largest absolute
%   value of the current x axis.
%
%   See also CSYM, YSYM.

if nargin < 1, m = max(abs(xlim)); end
xlim([-m,m])
