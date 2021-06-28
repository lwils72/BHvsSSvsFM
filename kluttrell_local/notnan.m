function [fnonans] = notnan(f)
% [fnonans] = notnan(f)
%   input: f is a column vector, or cell array of
%     column vectors, which may contain nans
%   output: fnonans is f with the nan values omitted 
%     (column vectors are different length than input)
%


  if iscell(f)==1
    inan=cellfun(@find,(cellfun(@isnan,f,'UniformOutput',false)),'UniformOutput',false);
    fnonans=f;
    for k=1:numel(f)
      fnonans{k}(inan{k})=[];
    end
  else
    fnonans=f;
    fnonans(find(isnan(f)))=[];
  end
