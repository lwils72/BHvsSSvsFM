function  [ff]=bandpass(f,dt,Tc0d,Tc1m,n)

%
% function  [ff]=bandpass(f,dt,Tc0d,Tc1m)
% wrapper to bandpass filter a dataset using
% a second order butterworth filter with given frequencies
%
%  f = time series to filter
%  dt = samplesize of f in seconds
%  Tc0d = highpass cuttoff period in days
%  Tc1m = lowpass cutoff period in minutes
%  [n] = [optional] order of filter (steepness of cutoff)
%

  if nargin<5,n=2;end
%
% design a highpass Butterworth filter
%
  flo=1/(Tc0d*24*3600); %lowfreq cutoff in Hz
  fhi=1/(Tc1m*60);  %highfreq cutoff in Hz
  fc=[flo fhi]; % cutoff frequency in Hz
  Wn=fc*2*dt; % normalized cutoff frequency (Wn=1 would be the nyquist, 2/dt Hz)
  %n=2; % n-th order filter (steepness of cutoff)

  [b,a]=butter(n,Wn);
%
% apply the zero-phase filter
%
  ff=filtfilt(b,a,f);
