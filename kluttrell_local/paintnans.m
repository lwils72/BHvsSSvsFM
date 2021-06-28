function [nf]=paintnans(f)
%
% [nf]=paintnans(f)
% interpolate the nans of the function f and return:
%   nf - the interpolated dataset with no nans
%
% This is similar to interp1 linear, but much faster.


%
  
  i=find(isnan(f));           % these are all the nan points
  nf=f; %this is the variable that will be transformed to nanless

  if numel(i)==0, return,end % in case there aren't any nans, just return f
  if numel(i)==numel(f), return,end % in case it's all nans, just return f

%
% fill in single points
% 
  if i(1)==1;i(1)=[];end % just in case the first point is nan 
  if numel(i)~=0, % in case only the first point was a NaN
    if i(numel(i))==numel(f);i=i(1:numel(i)-1);end % just in case the last point is nan 
    nf(i)=(f(i-1)+f(i+1))/2; %fill in the single points as neighbor average
  end
%
% fill in longer gaps
%  
  i=find(isnan(nf)); %these are now the nan segments
  if numel(i)==0, return,end % in case there aren't any nans, just return at this point
  
  di=diff(i);
  iskip=find(di>1);
  
  i1=i(iskip); % location of last nan (except last)
  i0=i(iskip+1); % location of first nan (except first)
  
  i0=[i(1);i0]; %location of all first nans
  i1=[i1;i(numel(i))]; % location of all last nans
  nskip=i1-i0+1; % number of nans per gap
  
  for j=1:numel(i0)
    if i0(j)==1;% just in case the first point is nan
      nf(i0(j):i1(j))=f(i1(j)+1);
    elseif i1(j)==numel(f);% just in case the last point is nan
      nf(i0(j):i1(j))=f(i0(j)-1);
    else
      nf(i0(j):i1(j))=linspace(f(i0(j)-1),f(i1(j)+1),nskip(j));
    end
  end
  
%
% only problem with this is that after a long gap
% sometimes the data is bad for a few samples... this doesn't get
% rid of the bad data... but it's good enough for now
%
