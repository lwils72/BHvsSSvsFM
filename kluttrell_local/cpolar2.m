function c = cpolar2(m)
%CPOLAR2   Wrapped version of cpolar: Shades of black, blue, white, and red color map
%   CPOLAR2(M) returns an M-by-3 matrix containing a "polar" colormap.
%   CPOLAR2, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   NOTE: m will be forced to be odd (4n+1) so that center value is white
%   Use with CSYM to ensure center value is zero.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(cpolar2)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT, CSYM, CPOLAR.


if nargin < 1, m = size(get(gcf,'colormap'),1); end

if mod(m,4)~=1, m=ceil(m/4)*4+1; end

  CR=[zeros(1,(m-1)/4),linspace(0,1,(m-1)/4+1),ones(1,(m-1)/4-1),linspace(1,0,(m-1)/4+1)]';
  
  CB=flipud(CR);

  CG=[zeros(1,(m-1)/4),linspace(0,1,(m-1)/4+1)]';
  CG=[CG;flipud(CG(1:end-1))];
  
  c=[CR,CG,CB];
