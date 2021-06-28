function y = ysym(m)
%YSYM   Symmetric y-axis about zero
%   YSYM(M) changes y axis limits to [-M,M].
%   YSYM, by itself, chooses M to be the largest absolute
%   value of the current y axis.
%
%   See also CSYM, XSYM.

if nargin < 1, m = max(abs(ylim)); end
ylim([-m,m])
