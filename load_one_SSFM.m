clear
addpath kluttrell_local
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

%
%   Load L&P station results and L&P EQ station pair Fast
%   Directions
%
    BHinfile='dats/Table_S1.csv';
    FMinfile='dats/FM_subsets_asof_20210615.shmax';
    SSinfile='dats/LiPeng2017tableS3.csv';
%   [B,S,L,F,Regions,P]=load_and_parse_SSBHFM(BHinfile,FMinfile,SSinfile);
%   unsure if this is needed
    
%
%   This excludes SS that are not in the LA basin 
%
     iLA=find(L.Y>P.R2(3) & L.Y<P.R2(4) & L.X>P.R2(1) & L.X<P.R2(2)); %list of iLA
     L.stationlon=L.stationlon(iLA);
     L.stationlat=L.stationlat(iLA);
     L.station=L.station(iLA);
     L.Nmeasurements=L.Nmeasurements(iLA);
     L.ResultantLength=L.ResultantLength(iLA);
     L.FastDirection=L.FastDirection(iLA);
     L.DelayTime=L.DelayTime(iLA);
     L.DelayTimeSTD=L.DelayTimeSTD(iLA);
     L.X=L.X(iLA);
     L.Y=L.Y(iLA);
     S.L_SHmax=S.L_SHmax(iLA,:,:);
     S.L_dSHmax=S.L_dSHmax(iLA,:,:);

%
%   Pick a specific station and try to reproduce their mean fast directions
%   so identify all EQs within the shear wave window under seismometer 
%   REMEMBER- Shear wave windows: epicenter distance less than hypocenter depth.
%

