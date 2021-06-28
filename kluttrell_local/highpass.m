function  [ff]=highpass(f,dt,Tcd,n)

%
% function  [ff]=highpass(f,dt,Tcd,n)
% wrapper to highpass filter a dataset using
% a first order butterworth filter with given frequency
%
%  f = time series to filter
%  dt = samplesize of f in seconds
%  Tcd = cuttoff period in days
%  [n] = [optional] order of filter (steepness of cutoff)
%

  if nargin<4,n=1;end
%
% design a highpass Butterworth filter
%
  fc=1/(Tcd*24*3600); % cutoff frquency in Hz
  Wn=fc*2*dt; % normalized cutoff frequency (Wn=1 would be the nyquist, 2/dt Hz)
  %n=1; % n-th order filter (steepness of cutoff)

  [b,a]=butter(n,Wn,'high');
%
% apply the zero-phase filter
%
  ff=filtfilt(b,a,f);
