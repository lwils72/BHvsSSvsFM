function Fpeak=peakdetect(f,P,Nsmooth,iqrsc)
% Fpeak=peakdetect(f,P,Nsmooth,iqrsc)
%   f = input frequency vector
%   P = input periodogram vector
%   Nsmooth = optional number of smoothing points 
%             for first derivative (default 10)
%   iqrsc = optional scaling of second derivative iqr
%           for the threshold of peak detection (default 2)
%   Fpeak = ouput frequency of detected peaks,
%           NaN if no peaks detected
%

if nargin < 4, Nsmooth=10;iqrsc=2;end

%
% detrend 
%
  dP=detrend(P);
%
% first derivative
%
  d1P=diff(dP);
%
% smoothed first derivative
%
  [b,a]=butter(1,2/Nsmooth,'low');
  d1Ps=filtfilt(b,a,d1P);
%
% (second) derivative of smoothed first derivative
%
  d2Ps=diff(d1Ps);
%
% find second derivatives below threshold (iqr*scale)
%
  i=find(d2Ps<-iqr(d2Ps)*iqrsc);
%
% identify freqs that belong to individual peaks
% ignore a "peak" with only one point below threshold
%
  R=contiguous(diff(i),1);R=R{2};
  if numel(R)==0,
    Fpeak=[NaN];
  else
    R(:,2)=R(:,2)+1;R=i(R)+1;
    if size(R,2)==1,R=R';end
%
% Identify the peak value within each peak group
%
    Fpeak=[];
    for k=1:size(R,1)
      [junk,j]=max(P(R(k,1):R(k,2)));
      Fpeak=[Fpeak;f(R(k,1)+j-1)];
    end

  end
