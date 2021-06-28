function ranout=ranany(binvalues,pdfvalues,N)

% ranout=ranany(binvalues,pdfvalues,N)
%
% Samples drawn from any distribution, specified by you.
% inputs:
%  - binvalues = row vector of centers of bins, usually evenly spaced but don't have to be, 1xM
%  - pdfvalues = probability corresponding to each bin, same size as binvalues, 1xM.
%      usually this sums to 1, but ranany will normalize just in case
%  - N = how many samples would you like to draw?
% outputs:
%  - ranout = Nx1 column vector of random numbers drawn from your distribution
%



pdfvalues=pdfvalues./sum(pdfvalues); % Normalize the PDF, if it's not already
cdfvalues=cumsum(pdfvalues,2); % calculate the CDF

rcdf=rand(N,1);
% [~,ibin]=min(abs(cdfvalues-rcdf),[],2); % hmmm, this one picks the "closest" bin value, but the resulting pdf is shifted half a bin negative
ibin=sum(cdfvalues<rcdf,2)+1; % this one picks the bin "above" the randomly selected cdf value, and seems to fix the problem
ranout=binvalues(ibin)';

