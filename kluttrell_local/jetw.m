function C = jetw(m,q)
%JETW   same as jet, but with exact middle shade white
%   JETW(M) returns an M-by-3 matrix containing a "polar" jet colormap.
%   JETW, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   JETW(M,Q) will set the center color to the gray shade defined by Q.
%   By default, Q=1 (white). Should be between [0,1].
%
%   NOTE: m will be forced to be odd to that center value is white
%   Use with CSYM to ensure center value is zero.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(jetw)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT, CSYM, CPOLAR.

%
% Version history
% 24 April 2012: originally created by Karen Luttrell, based on jet
%   - center color white
% 10 June 2015: revised to accept a second argument to allow
%     the center color to be defined (white default)
%


if nargin < 1, m = size(get(gcf,'colormap'),1); end
if nargin < 2, q = 1; end

if mod(m,2)==0, m=m+1; end
 
  m2=(m+1)/2;
  

  C=jet(m);
  
  C(m2,:)=[1 1 1]*q; %neutral gray shade for the exact middle
  
  
