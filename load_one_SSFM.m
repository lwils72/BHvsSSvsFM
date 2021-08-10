clear
addpath kluttrell_local
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
tic % start the clock to see how long the script takes to run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %   Load L&P station results and FM inversion results (and BH and other info)
  %     - L&P EQ station pair Fast Directions are loaded below
  %
    BHinfile='dats/Table_S1.csv';
    FMinfile='dats/FM_subsets_asof_20210615.shmax';
    SSinfile='dats/LiPeng2017tableS3.csv';
    [B,S,L,F,Regions,P]=load_and_parse_SSBHFM(BHinfile,FMinfile,SSinfile);
    
  %
  %   This excludes Seismic Stations that are not in the LA basin 
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
  % Load and look at whatever is available From the Li and Peng supplemnt,
  % These are the eq-station pairs we will use to try to replicate station averages
  %
    FORMAT=['%s %s ',repmat('%f ',1,23),repmat('%s ',1,4),'%f %s ',repmat('%f ',1,4),'%s ',repmat('%f ',1,5)];

    % This version loads only the high-qulaity eq-station pairs
    FILENAME='dats/LiPeng2017tableS2.csv'; %HiQuality SWS measurements (90k of them)
    fid=fopen(FILENAME);
    C=textscan(fid,FORMAT,'headerlines',1,'delimiter',',');
    fclose(fid);

    Leq.eventid=C{1}; % e.g., 15199577.CE.11369.HN
    Leq.station=C{2}; % e.g., 11369
    Leq.stationlat=C{3}; % e.g., 33.037
    Leq.stationlon=C{4}; % e.g., -115.624
    Leq.cuspid=C{5}; % e.g., 15199577
    Leq.year=C{6}; % e.g., 2012
    Leq.doy_det=C{7}; % e.g., 239.803
    Leq.eventlat=C{8}; % e.g., 33.0091
    Leq.eventlon=C{9}; % e.g., -115.553
    Leq.dist_event2station=C{10}; % e.g., 7.27946
    Leq.eventdepth=C{11}; % e.g., 9.064
    Leq.eventmag=C{12}; % e.g., 3.95
    Leq.backaz=C{13}; % e.g., 115.173
    Leq.spol=C{14}; % e.g., 96.665
    Leq.Dspol=C{15}; % e.g., 6.090
    Leq.wbeg=C{16}; % e.g., 4.414999
    Leq.wend=C{17}; % e.g., 5.704518
    Leq.dist_ruap_km=C{18}; % e.g., 
    Leq.dist_ruap_deg=C{19}; % e.g., 
    Leq.SNR=C{20}; % e.g., 21.0769
    Leq.tlag=C{21}; % e.g., 0.282500
    Leq.Dtlag=C{22}; % e.g., 0.005625
    Leq.fast=C{23}; % e.g., -21
    Leq.Dfast=C{24}; % e.g., 7.250
    Leq.anginc=C{25}; % e.g., 37.8 
    Leq.anginc_corr=C{26}; % e.g., anginc_corr
    Leq.type_ini=C{27}; % e.g., ass_5_16
    Leq.time=C{28}; % e.g., Sep
    Leq.comment=C{29}; % e.g., comment
    Leq.nyquist=C{30}; % e.g., 50
    Leq.gradeABCNR=C{31}; % e.g., ACl
    Leq.filt_lo=C{32}; % e.g., 
    Leq.filt_HI=C{33}; % e.g., 
    Leq.spolfastdiff=C{34}; % e.g., 62.335
    Leq.bandang=C{35}; % e.g., -29.6008
    Leq.pickgrade=C{36}; % e.g., UNDEFINE
    Leq.lambdamax=C{37}; % e.g., 3.6068532
    Leq.ndf=C{38}; % e.g., 25
    Leq.lambda2_min=C{39}; % e.g., 0.2680662E+06
    Leq.ttime=C{40}; % e.g., 3.69
    Leq.maxfreq=C{41}; % e.g., 4.45822

    % This version loads all available eq-station pairs, regardless of quality  
    FILENAME='dats/LiPeng2017tableS1.csv'; %RAW SWS measurements (230k of them)
    [~,r]=system(['cat ',FILENAME,' | grep -v cant']);
    C=textscan(r,FORMAT,'headerlines',1,'delimiter',',');

    Leq_all.eventid=C{1}; % e.g., 15199577.CE.11369.HN
    Leq_all.station=C{2}; % e.g., 11369
    Leq_all.stationlat=C{3}; % e.g., 33.037
    Leq_all.stationlon=C{4}; % e.g., -115.624
    Leq_all.cuspid=C{5}; % e.g., 15199577
    Leq_all.year=C{6}; % e.g., 2012
    Leq_all.doy_det=C{7}; % e.g., 239.803
    Leq_all.eventlat=C{8}; % e.g., 33.0091
    Leq_all.eventlon=C{9}; % e.g., -115.553
    Leq_all.dist_event2station=C{10}; % e.g., 7.27946
    Leq_all.eventdepth=C{11}; % e.g., 9.064
    Leq_all.eventmag=C{12}; % e.g., 3.95
    Leq_all.backaz=C{13}; % e.g., 115.173
    Leq_all.spol=C{14}; % e.g., 96.665
    Leq_all.Dspol=C{15}; % e.g., 6.090
    Leq_all.wbeg=C{16}; % e.g., 4.414999
    Leq_all.wend=C{17}; % e.g., 5.704518
    Leq_all.dist_ruap_km=C{18}; % e.g., 
    Leq_all.dist_ruap_deg=C{19}; % e.g., 
    Leq_all.SNR=C{20}; % e.g., 21.0769
    Leq_all.tlag=C{21}; % e.g., 0.282500
    Leq_all.Dtlag=C{22}; % e.g., 0.005625
    Leq_all.fast=C{23}; % e.g., -21
    Leq_all.Dfast=C{24}; % e.g., 7.250
    Leq_all.anginc=C{25}; % e.g., 37.8 
    Leq_all.anginc_corr=C{26}; % e.g., anginc_corr
    Leq_all.type_ini=C{27}; % e.g., ass_5_16
    Leq_all.time=C{28}; % e.g., Sep
    Leq_all.comment=C{29}; % e.g., comment
    Leq_all.nyquist=C{30}; % e.g., 50
    Leq_all.gradeABCNR=C{31}; % e.g., ACl
    Leq_all.filt_lo=C{32}; % e.g., 
    Leq_all.filt_HI=C{33}; % e.g., 
    Leq_all.spolfastdiff=C{34}; % e.g., 62.335
    Leq_all.bandang=C{35}; % e.g., -29.6008
    Leq_all.pickgrade=C{36}; % e.g., UNDEFINE
    Leq_all.lambdamax=C{37}; % e.g., 3.6068532
    Leq_all.ndf=C{38}; % e.g., 25
    Leq_all.lambda2_min=C{39}; % e.g., 0.2680662E+06
    Leq_all.ttime=C{40}; % e.g., 3.69
    Leq_all.maxfreq=C{41}; % e.g., 4.45822

  %
  % Older version: to get them all loaded quickly, just pull out what we want
  %

    % FILENAME='dats/LiPeng2017tableS1.csv'; %RAW SWS measurements (230k of them)
    % [~,r]=system(['cat ',FILENAME,' | grep -v cant | sed ''1d'' | awk -F, ''{print $4,$3,$9,$8,$11,$23}''']);
    % C=textscan(r,'%f %f %f %f %f %f');
    % 
    % Leq.stationlon=C{1}; % e.g., -115.624
    % Leq.stationlat=C{2}; % e.g., 33.037
    % Leq.eventlon=C{3}; % e.g., -115.553
    % Leq.eventlat=C{4}; % e.g., 33.0091
    % Leq.eventdepth=C{5}; % e.g., 9.064
    % Leq.fast=C{6}; % e.g., -21
       
  %
  % Recall: we alrady loaded the per-station "summary" fast directions above, as L
  %  - for the original proposal figure:
  %    Track the Y&H FM SHmax directions, at the location of the seismic stations
  %  - Only, this pulls out all of them, so also select only the ones from within the LA area
  %
    [~,r]=system('cat dats/LiPeng2017tableS3.csv | sed ''1d'' | awk -F, ''{print $2,$3}'' | gmt grdtrack -Ggrds/Yang.FMSHmax.grd -Z');
    C=textscan(r,'%f');
    L.FMSHmax=C{1};
    L.FMSHmax=L.FMSHmax(iLA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CAN WE REPRODUCE L&P's PER-STATION AVERAGE FAST DIRECTIONS FROM THE INDIVIDUAL EQ-STATION PAIRS?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  %   Pick a specific station and try to reproduce their mean fast directions
  %   so identify all EQs within the shear wave window under seismometer 
  %   REMEMBER- Shear wave windows: epicenter distance less than hypocenter depth.
  %

  %
  % Plot the station locations, along with their "number" so we can pick out a representative one
  %
    figure(1),clf,plot3(L.stationlon,L.stationlat,1:numel(L.station),'.k')
    view(2) % this sets the 3-D rotation view to be from above, the "x-y view"
    hold on,plot(P.c(:,1),P.c(:,2),'k') % plot the coast line
    axis(P.R3) % set the axis limits to focus on the LA area
    set(gca,'dataaspectratio',[1/cosd(mean(L.stationlat)) 1 1]) % set the relative size of the "x" and "y" units, so the map doesn't look squashed

  %
  % pick a single station
  %  - which "number" is it in our list?
  %  - what's is official station "name"?
  %  - what is it's estimated mean Fast Direction and uncertainty (from L&P, stored in L.)?
  %  - which earthquakes were recorded at this station that went into the 
  %
    figure(98),clf,plot(1:numel(L.X),L.Nmeasurements,'.-'),grid
    xlabel('station number (LA region only)'),ylabel('Nmeasurements: should be # eqs at each station (but isnt)')

    iOneSS=16; % to change which station we're looking at
     for k=1:numel(L.X) % use this code if you want to loop over all stations and look at the map one at a time
       iOneSS=k;        % use this if you want to loop

    OneSSname=L.station{iOneSS}; % JNH2 code for station 177 (from the whole). PASA for station 86 (from just LA ones). 
    OneSSfast=L.FastDirection(iOneSS); % the fast direction

    ieqOneSS=find(strcmp(Leq.station,OneSSname)); % which eqs (from the high-quality list) were recorded at this station?
    Leq.fast(ieqOneSS); % fast direction of those 20 eqs (not stored here, just for our reference)
    % Note: already I can tell there's something wrong.  According to LiPengTableS3, the agregate fast direction for 
    %   station PASA has 18 measurements, but from LiPengTableS2 (the high quality ones), it only has one eq listed.
    %   Try the longer list from LiPengTableS3, see if there are 18 PASA eqs in that list...

    ieqOneSS_all=find(strcmp(Leq_all.station,OneSSname)); % which eqs (from the all-quality list) were recorded at this station?
    Leq_all.fast(ieqOneSS); % fast direction of those 20 eqs (not stored here, just for our reference)
    % Nope, this has 7 eqs, some are duplicated but still have different fast directions indicated...
    %  - so in what sense does PASA have 18 "measurements"?  No idea, and text gives no indication
    %  - Try a different station other than PASA... does it add up? Nope, L.Nmeasurements way over-reports the number of eq-station pairs
    %    (no indication in tables or text that each eq corresponds to multiple measurements, e.g.)
    %  - OK, after trying several of these, the formula seems to be usually Nmeasurements = 2*Neq, and sometimes N > 2*Neq... 
    %    mysterious... but now we can move on, assuming "Nmeasurements" is unreliable...
 
    ieqOneSS_shallow=find(strcmp(Leq_all.station,OneSSname)&(Leq.eventdepth<=5))%list of SS with only shallow eq
    ieqOneSS_deep=find(strcmp(Leq_all.station,OneSSname)&(Leq.eventdepth>5))%list of SS with eq that are larger than 5km
    
    %stop
    
    [Leq_all.eventlon(ieqOneSS_all),Leq_all.eventlat(ieqOneSS_all),Leq_all.eventdepth(ieqOneSS_all),Leq_all.eventmag(ieqOneSS_all),...
     Leq_all.fast(ieqOneSS_all),Leq_all.Dfast(ieqOneSS_all)]; %prints the eqs out so we can see them... when not semi-coloned

  %
  % Plot these earthquakes on the map, and the station they're recorded at
  %
    figure(1),
    plot(L.stationlon(iOneSS),L.stationlat(iOneSS),'go','linewidth',1)
    plot(Leq.eventlon(ieqOneSS),Leq.eventlat(ieqOneSS),'r*')
    plot(Leq_all.eventlon(ieqOneSS_all),Leq_all.eventlat(ieqOneSS_all),'bs')
    title(['station ',OneSSname,' and ',num2str(L.Nmeasurements(iOneSS)),' measurements from ',num2str(numel(ieqOneSS_all)),' earthquakes'])
    legend('stations','coast','One Station','HighQ Eqs','All Eqs','location','southwest')

  %
  % plot a histogram of the fast direction population at our One Station,
  % Compare with the Li&Peng mean fast direction.
  % Do these each make sense? (flip through the stations and check)
  %  - ANSWER: mostly yes.  It seems that these fast directions are actually calculated using
  %    the full list of eqs, not just the high quality ones (as was stated in the text).
  %    See, e.g., station #17 ("24853"), where the mean direction is clearly the mean of "all", not "high Q"
  %    (station #21 is a good example of very poorly constrained fast direction...)
  %
    binedges=[-90:10:90]-1e-4; % shifting bins a smidgen makes it behave better for angles right on the edge of a bin (e.g., 30deg)

    figure(2),clf
    subplot(211)
      histogram(Leq_all.fast(ieqOneSS_all),binedges)
      hold on
      histogram(Leq.fast(ieqOneSS),binedges)
      xlabel('Fast Direction Degrees')
      ylabel('Number of EQs pairs')
      title(['Histogram Fast Directions for station ',OneSSname])
      plot(L.FastDirection(iOneSS)*[1 1],ylim,'k','linewidth',2)
      legend('all eqs','high Q eqs','LiPeng meanfast','location','northwest')
      xlim([-90,90])
      xticks(-90:15:90)
      grid
    subplot(223)
      polarhistogram(deg2rad(Leq_all.fast(ieqOneSS_all)),deg2rad(binedges),'facecolor','k');hold on
      polarhistogram(deg2rad(Leq_all.fast(ieqOneSS_all))+pi,deg2rad(binedges)+pi,'facecolor','k'); % plot the data twice, for 180º symmetry
      polarplot([1 1]*deg2rad(L.FastDirection(iOneSS)),[-1 1]*max(rlim),'r','linewidth',1)
      title('all eqs')
      set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise') % make axes degEofN, not degNofE
    subplot(224)
      polarhistogram(deg2rad(Leq.fast(ieqOneSS)),deg2rad(binedges),'facecolor','k');hold on
      polarhistogram(deg2rad(Leq.fast(ieqOneSS))+pi,deg2rad(binedges)+pi,'facecolor','k'); % plot the data twice, for 180º symmetry
      polarplot([1 1]*deg2rad(L.FastDirection(iOneSS)),[-1 1]*max(rlim),'r','linewidth',1)
      title('high quality eqs')
      set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise') % make axes degEofN, not degNofE
      
     figure(21),clf
     subplot(211)
      histogram(Leq_all.fast(ieqOneSS_shallow),binedges)
      hold on
      histogram(Leq.fast(ieqOneSS),binedges)
      xlabel('Fast Direction Degrees')
      ylabel('Number of SHALLOW EQs pairs')
      title(['Histogram Fast Directions for station ',OneSSname])
      plot(L.FastDirection(iOneSS)*[1 1],ylim,'k','linewidth',2)
      legend('all eqs','high Q eqs','LiPeng meanfast','location','northwest')
      xlim([-90,90])
      xticks(-90:15:90)
      grid
    subplot(223)
      polarhistogram(deg2rad(Leq_all.fast(ieqOneSS_shallow)),deg2rad(binedges),'facecolor','k');hold on
      polarhistogram(deg2rad(Leq_all.fast(ieqOneSS_shallow))+pi,deg2rad(binedges)+pi,'facecolor','k'); % plot the data twice, for 180º symmetry
      polarplot([1 1]*deg2rad(L.FastDirection(iOneSS)),[-1 1]*max(rlim),'r','linewidth',1)
      title('all eqs')
      set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise') % make axes degEofN, not degNofE
    subplot(224)
      polarhistogram(deg2rad(Leq.fast(ieqOneSS)),deg2rad(binedges),'facecolor','k');hold on
      polarhistogram(deg2rad(Leq.fast(ieqOneSS))+pi,deg2rad(binedges)+pi,'facecolor','k'); % plot the data twice, for 180º symmetry
      polarplot([1 1]*deg2rad(L.FastDirection(iOneSS)),[-1 1]*max(rlim),'r','linewidth',1)
      title('high quality eqs')
      set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise') % make axes degEofN, not degNofE
      
      figure(22),clf
    subplot(211)
      histogram(Leq_all.fast(ieqOneSS_deep),binedges)
      hold on
      histogram(Leq.fast(ieqOneSS),binedges)
      xlabel('Fast Direction Degrees')
      ylabel('Number of DEEP EQs pairs')
      title(['Histogram Fast Directions for station ',OneSSname])
      plot(L.FastDirection(iOneSS)*[1 1],ylim,'k','linewidth',2)
      legend('all eqs','high Q eqs','LiPeng meanfast','location','northwest')
      xlim([-90,90])
      xticks(-90:15:90)
      grid
    subplot(223)
      polarhistogram(deg2rad(Leq_all.fast(ieqOneSS_deep)),deg2rad(binedges),'facecolor','k');hold on
      polarhistogram(deg2rad(Leq_all.fast(ieqOneSS_deep))+pi,deg2rad(binedges)+pi,'facecolor','k'); % plot the data twice, for 180º symmetry
      polarplot([1 1]*deg2rad(L.FastDirection(iOneSS)),[-1 1]*max(rlim),'r','linewidth',1)
      title('all eqs')
      set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise') % make axes degEofN, not degNofE
    subplot(224)
      polarhistogram(deg2rad(Leq.fast(ieqOneSS)),deg2rad(binedges),'facecolor','k');hold on
      polarhistogram(deg2rad(Leq.fast(ieqOneSS))+pi,deg2rad(binedges)+pi,'facecolor','k'); % plot the data twice, for 180º symmetry
      polarplot([1 1]*deg2rad(L.FastDirection(iOneSS)),[-1 1]*max(rlim),'r','linewidth',1)
      title('high quality eqs')
      set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise') % make axes degEofN, not degNofE

  %
  % Calculate the mean fast direction from scratch, from the "all eqs" list of fast directions
  % Compare to the L&P calculated directions and population... do they all make sense?
  %  (station #16 is a good test case and demo of "regular mean gives the wrong answer", or #52 is even simpler with just 2 angles)
  %  - ANSWER: The new functions defined below work pretty well, very close to L&P answers most of the time, but
  %    sometimes really different (e.g., station # 16... sigh... still doesn't work well...)
  %    (maybe the difference is in we're neglecting the weights, but sometimes maybe it's not...)
  %    E.g., for station #38, their average is clearly not sensible...
  % 
  %
    % define new funtions that will take the 360º or 180º wrapped circular
    % mean of directional data, one line subroutine to repeat (inline
    % define funtion)
    mean360= @(x) atan2d(mean(sind(x)),mean(cosd(x))); % our very own regular circular mean...
    mean180= @(x) atan2d(mean(sind(x*2)),mean(cosd(x*2)))/2; % adapted for 180º wrapping
    mean180weighted= @(x,y) atan2d(sum(sind(x.*y*2))/sum(y),sum(cosd(x.*y*2))/sum(y))/2; % adapted for 180º wrapping, and weights

    OurMeanFast=mean180(Leq.fast(ieqOneSS));
    OurMeanFast_all=mean180(Leq_all.fast(ieqOneSS_all));

    OurMeanFast_all_weighted=mean180weighted(Leq_all.fast(ieqOneSS_all),Leq_all.tlag(ieqOneSS_all)); % maybe weight by time lag?

    subplot(223)
      polarplot([1 1]*deg2rad(OurMeanFast_all),[-1 1]*max(rlim),'b','linewidth',1)
    subplot(224)
      polarplot([1 1]*deg2rad(OurMeanFast),[-1 1]*max(rlim),'b','linewidth',1)


     pause % use these to end the loop started above, if you want to look at each station one at a time
     end

  %
  % OK, AT THIS POINT, we're pretty sure we're calculating the 180º circular mean correctly, but not getting their answer.
  % We can either continue to try to figure out what they did (and maybe email them and ask), or we can move on and
  % continue testing if we get different answers for fast direction from eqs from different depths... (probably the latter)
  %



toc % stop the clock and report how long it took
stop

%
% OLDER STUFF
%  
  %
  %   to use circ_mean need to convert angles to radians
  %

  % r2d(angle(mean(exp(d2r(L.FastDirection*2)*i))))/2
  % sum(cosd(L.FastDirection))=y
  % sum(sind(L.FastDirection))=x
  % THEmean=atand(sum(sind(L.FastDirection*2)/sum(cosd(L.FastDirection*2)))) %atand(x/y)  

  %
  %   Possibly generate histogram to show the 20 fast directions vs the mean
  %   plot various results for each of the pairs L(summary nums) vs Leq(all data)
  %   play around see what all fits and makes sense... 
  %   calculate the circular mean fast direction 
  %


  %
  %   Only graphing the fast directions and the single mean 
  %

  % figure(61),clf
  % plot(L.FastDirection(iOneSS),'^','linewidth',1)
  % plot(Leq.fast,'^','linewidth',1)
  % axis equal,xlim([-30,120]),ylim([-90,90]),grid
  % xticks((-30:30:120)),yticks((-90:30:90))
  % 
  % figure(62),clf
  % plot(L.FastDirection(iOneSS),'^','linewidth',1)
  % plot(Leq.fast(iOneSS),'^','linewidth',1)
  % axis equal,xlim([-30,120]),ylim([-90,90]),grid
  % xticks((-30:30:120)),yticks((-90:30:90))
  % 
  %  figure(70),clf
  %  plot(L.FastDirection,Leq.fast(iOneSS),'^','linewidth',1)
  %   
  %  figure(71),clf
  %  plot(L.FastDirection,Leq.fast(iOneSS),'^','linewidth',1)
  % 
  % figure(80),clf
  % plot(Leq.stationlon,Leq.stationlat,'k^')
  % hold on, plot(L.X(:,1),L.Y(:,2),'-k')
  % scatter(Leq.eventlon,L.eventlat,5,Leq.eventdepth,'filled'),colorbar

