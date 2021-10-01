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
   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %     Comparing Mean Fast Directions from Li and Peng vs our Calculated 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    
     for k=1:numel(L.X) % use this code if you want to loop over all stations and look at the map one at a time
       iOneSS=k;        % use this if you want to loop
       
       
    OneSSname=L.station{iOneSS}; % JNH2 code for station 177 (from the whole). PASA for station 86 (from just LA ones). 
    OneSSfast=L.FastDirection(iOneSS); % the fast direction

    ieqOneSS=find(strcmp(Leq.station,OneSSname)); % which eqs (from the high-quality list) were recorded at this station?
    Leq.fast(ieqOneSS); % fast direction of those 20 eqs (not stored here, just for our reference)
    
    ieqOneSS_all=find(strcmp(Leq_all.station,OneSSname)); % which eqs (from the all-quality list) were recorded at this station?
    Leq_all.fast(ieqOneSS); % fast direction of those 20 eqs (not stored here, just for our reference)
       
    ieqOneSS_shallow=find(strcmp(Leq_all.station,OneSSname)&(Leq_all.eventdepth<=5)); %list of SS with only shallow eq
    ieqOneSS_deep=find(strcmp(Leq_all.station,OneSSname)&(Leq_all.eventdepth>5));%list of SS with eq that are larger than 5km
    
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
    MeanFastShallow(k,1)=mean180(Leq_all.fast(ieqOneSS_shallow));
    MeanFastDeep(k,1)=mean180(Leq_all.fast(ieqOneSS_deep));
    
     end
     
    

    %
    %   stations colored by the mean of the fast direction 
    %
    
    sc=10; %how many km long should the SHmax stick be?

    figure(1),clf
    %set(gcf,'Position',[100,400,665,440])
        subplot(3,3,1),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(MeanFastDeep),L.X+sc*sind(MeanFastDeep)]',[L.Y-sc*cosd(MeanFastDeep),L.Y+sc*cosd(MeanFastDeep)]','k','linewidth',2)
        scatter(L.X,L.Y,20,MeanFastDeep,'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'a) SWS Mean Fast Direction - Deep EQS','FontSize',12,'FontWeight','bold')
       
        subplot(3,3,2),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(MeanFastShallow),L.X+sc*sind(MeanFastShallow)]',[L.Y-sc*cosd(MeanFastShallow),L.Y+sc*cosd(MeanFastShallow)]','k','linewidth',2)
        scatter(L.X,L.Y,20,MeanFastShallow,'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'b) SWS Mean Fast Direction - Shallow EQS','FontSize',12,'FontWeight','bold')
        
        subplot(3,3,3),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',2)
        scatter(L.X,L.Y,20,L.FastDirection,'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'a) SS MeanFastAll','FontSize',12,'FontWeight','bold')
        
     %
     % stations colored by fmshmax fast directions 
     %
     
     idmax=find(S.dmaxset==6);
     izmaxshallow=find(S.zmaxset==5);
     izmaxdeep=find(S.zmaxset==999);
     
        subplot(3,3,4),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(S.SHmax(:,idmax,izmaxdeep)),L.X+sc*sind(S.SHmax(:,idmax,izmaxdeep))]',[L.Y-sc*cosd(S.SHmax(:,idmax,izmaxdeep)),L.Y+sc*cosd(S.SHmax(:,idmax,izmaxdeep))]','k','linewidth',2)
        scatter(L.X,L.Y,20,S.SHmax(:,idmax,izmaxdeep),'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'c) FM SHmax - All EQs','FontSize',12,'FontWeight','bold')

        subplot(3,3,5),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(S.SHmax(:,idmax,izmaxshallow)),L.X+sc*sind(S.SHmax(:,idmax,izmaxshallow))]',[L.Y-sc*cosd(S.SHmax(:,idmax,izmaxshallow)),L.Y+sc*cosd(S.SHmax(:,idmax,izmaxshallow))]','k','linewidth',2)
        scatter(L.X,L.Y,20,S.SHmax(:,idmax,izmaxshallow),'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'d) FM SHmax - Shallow EQs','FontSize',12,'FontWeight','bold')
        
        subplot(3,3,6),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(S.SHmax(:,idmax,izmaxdeep)),L.X+sc*sind(S.SHmax(:,idmax,izmaxdeep))]',[L.Y-sc*cosd(S.SHmax(:,idmax,izmaxdeep)),L.Y+sc*cosd(S.SHmax(:,idmax,izmaxdeep))]','k','linewidth',2)
        scatter(L.X,L.Y,20,S.SHmax(:,idmax,izmaxdeep),'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'b) FMSHmax All','FontSize',12,'FontWeight','bold')
        
     %   
     % stations colored by the circular mean difference 
     %
     
         %
         % 180 circular difference used to compare our fast directions within the
         % shallow and deep maps
         %
         
        daz = @(x,y) (mod(y-x+90,180)-90);  % Circular Difference
        rmsaz = @(x,y) sqrt(sum(daz(x,y).^2,'omitnan')./sum(~isnan(y))); % RMS of Circular Difference, ignores NaNs, works with arrays of column vectors
        %zBHvFM=abs(daz(S.SHmax,S.B_SHmax)); % Absolute Circular Difference of the Data
        FMFDDiffDeep=daz(MeanFastDeep,S.SHmax(:,idmax,izmaxdeep));
        FMFDDiffShallow=daz(MeanFastShallow,S.SHmax(:,idmax,izmaxshallow));
        FMFDDiffAll=daz(L.FastDirection,S.SHmax(:,idmax,izmaxdeep));
     
         subplot(3,3,7),
         plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
         hold on,
         scatter(L.X,L.Y,20,abs(FMFDDiffDeep),'filled'),colorbar
         xticks([]),yticks([])
         caxis([0,90])
         text(227,3737,'c) SS CD Deep','FontSize',12,'FontWeight','bold')

         subplot(3,3,8),
         plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
         hold on,
         scatter(L.X,L.Y,20,abs(FMFDDiffShallow),'filled'),colorbar
         xticks([]),yticks([])
         caxis([0,90])
         text(227,3737,'c) SS CD Shallow','FontSize',12,'FontWeight','bold')
         
         subplot(3,3,9),
         plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
         hold on,
         scatter(L.X,L.Y,20,abs(FMFDDiffAll),'filled'),colorbar
         xticks([]),yticks([])
         caxis([0,90])
         text(227,3737,'c) SS CD All','FontSize',12,'FontWeight','bold')
         stop
    %
    % 180 circular difference used to compare our fast directions within the
    % shallow and deep maps
    %
  
%stacklocalPDFs.m:    daz = @(x,y) (mod(y-x+90,180)-90);  % Circular Difference
%stacklocalPDFs.m:    rmsaz = @(x,y) sqrt(sum(daz(x,y).^2,'omitnan')./sum(~isnan(y))); % RMS of Circular Difference, ignores NaNs, works with arrays of column vectors
%stacklocalPDFs.m:    %zBHvFM=abs(daz(S.SHmax,S.B_SHmax)); % Absolute Circular Difference of the Data
%stacklocalPDFs.m:    L.SHmaxDiff=daz(L.SHmax1,L.SHmax2);
    
    
    %
    %   Plotting BH and SWS data to try to compare the heterogenious data 
    %
    
    figure(2),clf
        subplot(2,1,1)
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',2)
        scatter(L.X,L.Y,20,(L.FastDirection),'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'a) SWS FastDirections','FontSize',18,'FontWeight','bold')
        
        subplot(2,1,2)
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        hold on,
        plot([B.X-sc*sind(B.SHmax),B.X+sc*sind(B.SHmax)]',[B.Y-sc*cosd(B.SHmax),B.Y+sc*cosd(B.SHmax)]','k','linewidth',2)
        scatter(B.X,B.Y,20,(B.SHmax),'filled'),colorbar
        xticks([]),yticks([])
        caxis([-90,90])
        text(227,3737,'b) BH SHmax','FontSize',18,'FontWeight','bold')
        
        stop
       
       
    