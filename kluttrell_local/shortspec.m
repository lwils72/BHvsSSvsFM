function [F,T,P]=shortspec(g,t0,dt,winhr,foverlap,Fmax,Fmin)
%[F,T,P]=shortspec(g,t0,dt,winhr,foverlap,[Fmax],[Fmin])
%  wrapper for spectrogram subroutine
%
%  input:
%     g = time series
%     t0 = initial time (in MATLAB time)
%     dt = sampling interval (in seconds)
%     winhr = how many hours per window in time
%     foverlap = fraction of window overlap
%     [Fmax] = maximum frequency to return (in cycles per hour)
%     [Fmin] = minimum frequency to return (in cycles per hour)
%  output:
%     F = vector of frequencied (in cycles per hour)
%     T = vector of times (in MATLAB time)
%     P = log spectrogram, of F and T space
%


%
% e.g.
%
  %winhr=12;              % about how many hours per window in time
  %foverlap=7/8;         % fraction of window overlap
  %Fmax=2;    % 30 minutes
  %Fmin=60/300;  % 5 hours


%
%spectrogram parameters
%
  wind=winhr*3600/dt;
  noverlap=wind*foverlap;
  nfft=wind;

%
%spectrogram
%
  [S,F,T,P] = spectrogram(g,wind,noverlap,nfft,1/dt);
  
  
  P=log10(P);
  F=F*3600;
  T=T/3600/24+t0;
  
%
% trim to desired frequency range
%
  if nargin > 5
    P=P(find(F<=Fmax),:);
    F=F(find(F<=Fmax));
  end

  if nargin == 7
    P=P(find(F>=Fmin-eps),:);
    F=F(find(F>=Fmin-eps));
  end

