clear
addpath kluttrell_local 
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
%
% Load the earlier stuff.. actually there's alot of it...
%
  tic
  BHinfile='dats/Table_S1.csv';
  FMinfile='dats/FM_subsets_asof_20200108.shmax';
  [B,~,F,Regions,P]=load_and_parse_BHFM(BHinfile,FMinfile);


%
% Plot what we know
% Mapview plot of borehole locations, colored by number because it helps to easily identify them.
%
  figure(1),clf
  plot(P.cutm(:,1),P.cutm(:,2),'k'),hold on,
  plot(P.futm(:,1),P.futm(:,2),'color',[1 1 1]*0.75)
  axis equal,axis(P.R2)
  scatter3(B.X,B.Y,1:numel(B.X),50,1:numel(B.X),'filled'),colorbar
  plot3(B.X,B.Y,1:numel(B.X),'ok')
  colormap(jet),view(2) % sets view to X-Y 2D 
  set(gcf,'renderer','Painters')
  title('Boreholes, colored by number...')

%
% Load and look at whatever is available From the Li and Peng supplemnt
%
  % FILENAME='dats/LiPeng2017tableS1.csv'; %RAW SWS measurements (230k of them)
  % FORMAT=['%s %s ',repmat('%f ',1,23),repmat('%s ',1,4),'%f %s ',repmat('%f ',1,4),'%s ',repmat('%f ',1,5)];
  % [~,r]=system(['cat ',FILENAME,' | grep -v cant']);
  % C=textscan(r,FORMAT,'headerlines',1,'delimiter',',');

  % FILENAME='dats/LiPeng2017tableS2.csv'; %HiQuality SWS measurements (90k of them)
  % fid=fopen(FILENAME);
  % C=textscan(fid,FORMAT,'headerlines',1,'delimiter',',');
  % fclose(fid);

  % S.eventid=C{1}; % e.g., 15199577.CE.11369.HN
  % S.station=C{2}; % e.g., 11369
  % S.stationlat=C{3}; % e.g., 33.037
  % S.stationlon=C{4}; % e.g., -115.624
  % S.cuspid=C{5}; % e.g., 15199577
  % S.year=C{6}; % e.g., 2012
  % S.doy_det=C{7}; % e.g., 239.803
  % S.eventlat=C{8}; % e.g., 33.0091
  % S.eventlon=C{9}; % e.g., -115.553
  % S.dist_event2station=C{10}; % e.g., 7.27946
  % S.eventdepth=C{11}; % e.g., 9.064
  % S.eventmag=C{12}; % e.g., 3.95
  % S.backaz=C{13}; % e.g., 115.173
  % S.spol=C{14}; % e.g., 96.665
  % S.Dspol=C{15}; % e.g., 6.090
  % S.wbeg=C{16}; % e.g., 4.414999
  % S.wend=C{17}; % e.g., 5.704518
  % S.dist_ruap_km=C{18}; % e.g., 
  % S.dist_ruap_deg=C{19}; % e.g., 
  % S.SNR=C{20}; % e.g., 21.0769
  % S.tlag=C{21}; % e.g., 0.282500
  % S.Dtlag=C{22}; % e.g., 0.005625
  % S.fast=C{23}; % e.g., -21
  % S.Dfast=C{24}; % e.g., 7.250
  % S.anginc=C{25}; % e.g., 37.8 
  % S.anginc_corr=C{26}; % e.g., anginc_corr
  % S.type_ini=C{27}; % e.g., ass_5_16
  % S.time=C{28}; % e.g., Sep
  % S.comment=C{29}; % e.g., comment
  % S.nyquist=C{30}; % e.g., 50
  % S.gradeABCNR=C{31}; % e.g., ACl
  % S.filt_lo=C{32}; % e.g., 
  % S.filt_HI=C{33}; % e.g., 
  % S.spolfastdiff=C{34}; % e.g., 62.335
  % S.bandang=C{35}; % e.g., -29.6008
  % S.pickgrade=C{36}; % e.g., UNDEFINE
  % S.lambdamax=C{37}; % e.g., 3.6068532
  % S.ndf=C{38}; % e.g., 25
  % S.lambda2_min=C{39}; % e.g., 0.2680662E+06
  % S.ttime=C{40}; % e.g., 3.69
  % S.maxfreq=C{41}; % e.g., 4.45822

%
% OK, to get them all loaded, just pull out what we want
%

  FILENAME='dats/LiPeng2017tableS1.csv'; %RAW SWS measurements (230k of them)
  [~,r]=system(['cat ',FILENAME,' | grep -v cant | sed ''1d'' | awk -F, ''{print $4,$3,$9,$8,$11,$23}''']);
  C=textscan(r,'%f %f %f %f %f %f');

  S.stationlon=C{1}; % e.g., -115.624
  S.stationlat=C{2}; % e.g., 33.037
  S.eventlon=C{3}; % e.g., -115.553
  S.eventlat=C{4}; % e.g., 33.0091
  S.eventdepth=C{5}; % e.g., 9.064
  S.fast=C{6}; % e.g., -21

%
% Well, that's only reading in the first 1% of them, but go ahead and plot these anyways, to see what's in there
%

  figure(2),clf
  plot(S.stationlon,S.stationlat,'k^')
  hold on, plot(P.c(:,1),P.c(:,2),'-k')
  scatter(S.eventlon,S.eventlat,5,S.eventdepth,'filled'),colorbar
  axis([-121,-114,32,36]),grid
  set(gca,'dataaspectratio',[1/cosd(34),1, 1])
  title('seismic stations, earthquakes colored by depth (230k of them)')

  figure(3),clf
  plot(S.stationlon,S.stationlat,'k^')
  hold on, plot(P.c(:,1),P.c(:,2),'-k')
  scatter(S.eventlon,S.eventlat,5,S.fast,'filled'),colorbar
  axis([-121,-114,32,36]),grid
  set(gca,'dataaspectratio',[1/cosd(34),1, 1])
  colormap(hsv),caxis([-90,90])
  title('seismic stations, earthquakes colored by fast axis (230k of them)')

%
% Load the second version.  Actually, try just pulling out what we really want from this.  That loads them all...
%
  FILENAME='dats/LiPeng2017tableS2.csv'; %HiQuality SWS measurements (90k of them)
  [~,r]=system(['cat ',FILENAME,' | sed ''1d'' | awk -F, ''{print $4,$3,$9,$8,$11,$23}''']);
  C=textscan(r,'%f %f %f %f %f %f');

  T.stationlon=C{1}; % e.g., -115.624
  T.stationlat=C{2}; % e.g., 33.037
  T.eventlon=C{3}; % e.g., -115.553
  T.eventlat=C{4}; % e.g., 33.0091
  T.eventdepth=C{5}; % e.g., 9.064
  T.fast=C{6}; % e.g., -21

  [T.X,T.Y]=ll2utm(T.eventlat,T.eventlon);
  T.X=T.X/1e3;
  T.Y=T.Y/1e3;

  figure(12),clf
  plot(T.stationlon,T.stationlat,'k^')
  hold on, plot(P.c(:,1),P.c(:,2),'-k')
  scatter(T.eventlon,T.eventlat,5,T.eventdepth,'filled'),colorbar
  axis([-121,-114,32,36]),grid
  set(gca,'dataaspectratio',[1/cosd(34),1, 1])
  title('seismic stations, earthquakes colored by depth (90k of them)')

  figure(13),clf
  plot(T.stationlon,T.stationlat,'k^')
  hold on, plot(P.c(:,1),P.c(:,2),'-k')
  scatter(T.eventlon,T.eventlat,5,T.fast,'filled'),colorbar
  axis([-121,-114,32,36]),grid
  set(gca,'dataaspectratio',[1/cosd(34),1, 1])
  colormap(hsv),caxis([-90,90])
  title('seismic stations, earthquakes colored by fast axis (90k of them)')

%
% OK, now load the "summary" version
%
  fid=fopen('dats/LiPeng2017tableS3.csv');
  C=textscan(fid,'%s %f %f %f %f %f %f %f','headerlines',1,'delimiter',',');
  fclose(fid);

  L.station=C{1};
  L.stationlon=C{2};
  L.stationlat=C{3};
  L.Nmeasurements=C{4};
  L.ResultantLength=C{5};
  L.FastDirection=C{6}; % deg E of N
  L.DelayTime=C{7};
  L.DelayTimeSTD=C{8};

  % convert to UTM km
  [L.X,L.Y]=ll2utm(L.stationlat,L.stationlon);
  L.X=L.X/1e3;
  L.Y=L.Y/1e3;

  % Track the Y&H FM SHmax directions, at the location of the seismic stations
  [~,r]=system('cat dats/LiPeng2017tableS3.csv | sed ''1d'' | awk -F, ''{print $2,$3}'' | gmt grdtrack -Ggrds/Yang.FMSHmax.grd -Z');
  C=textscan(r,'%f');
  L.FMSHmax=C{1};



  figure(4),clf,
  plot(L.stationlon,L.stationlat,'k^')
  hold on, plot(P.c(:,1),P.c(:,2),'-k')
  scatter(L.stationlon,L.stationlat,50,L.FastDirection,'filled'),colorbar
  axis([-121,-114,32,36]),grid
  set(gca,'dataaspectratio',[1/cosd(34),1, 1])
  colormap(hsv),caxis([-90,90])
  title('seismic stations, colored by fast axis (90k of them)')

%
% OK, so let's start by just comparing where we have data...
%
  figure(5),clf
  plot(F.x,F.y,'.','color',[1 1 1]*0.75),hold on
  plot(B.x,B.y,'o',L.stationlon,L.stationlat,'^','linewidth',2)
  contour(P.CVM.x,P.CVM.y,P.CVM.z/1e3,(-20:2:1),'linewidth',1),colorbar,caxis([-20,1])
  plot(P.c(:,1),P.c(:,2),'-k')
  axis([-121,-114,32,36]),grid
  set(gca,'dataaspectratio',[1/cosd(34),1, 1])
  legend('FMs','BHs','SWSs')
  title('seismic stations, colored by fast axis (407 of them)')

  %
  % try making a depth color scheme that matches the one from GMT, based on magma
  %
    C=[0.9844    0.9883    0.7461;... % 1 km
       0.9922    0.7588    0.5278;... % 3 km
       0.9814    0.5293    0.3789;... % 5 km
       0.8975    0.3154    0.3892;... % 7 km
       0.7109    0.2129    0.4746;... % 9 km
       0.5088    0.1470    0.5039;... % 11 km
       0.3145    0.0703    0.4834;... % 13 km
       0.1128    0.0659    0.2759;... % 15 km
            0         0    0.0156;... % 17 km
            0         0    0.0156];   % 19 km (same as 17, just cause)
    Ninterp=20; % how many interpolations would you like between each?
    C=flipud(interp1(1:length(C),C,1:1/Ninterp:length(C)));
    C=[ones(Ninterp/2,1)*C(1,:);C;ones(Ninterp/2,1)*C(end,:)]; % add half a matching step on each end so it goes from 0 - 20 km
    colormap(C)

%
% OK, can I make a map with them colored by inferred SHmax, and by stick
% BH vs Y&H, SWS vs Y&H, map with BH and SWS, map with FM at those locations...
%  - that would be enough for the inital comparison, I think.
%
  sc=10; %how many km long should the SHmax stick be?
  
  figure(6),clf

  set(gcf,'Position',[100,351,877,489])
  [ha,pos]=tight_subplot(2,2,0.05,0.15,0.05); % this just helps it be a nicer looking journal figure
 
  axes(ha(1))
  % subplot(221)
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,
    plot([B.X-sc*sind(B.SHmax),B.X+sc*sind(B.SHmax)]',[B.Y-sc*cosd(B.SHmax),B.Y+sc*cosd(B.SHmax)]','k','linewidth',2)
    scatter(B.X,B.Y,50,B.SHmax,'filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'a) borehole SHmax','FontSize',12,'FontWeight','bold')

  axes(ha(2))
  % subplot(222)
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,

    scatter(T.X,T.Y,5,T.fast,'filled'),

    plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',2)
    scatter(L.X,L.Y,50,L.FastDirection,'^','filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'b) SWS fast direction','FontSize',12,'FontWeight','bold')
  
  axes(ha(3))
  % subplot(223)
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,
    plot([B.X-sc*sind(B.FMSHmax),B.X+sc*sind(B.FMSHmax)]',[B.Y-sc*cosd(B.FMSHmax),B.Y+sc*cosd(B.FMSHmax)]','k','linewidth',2)
    scatter(B.X,B.Y,50,B.FMSHmax,'filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'c) Y&H13 SHmax at BHs','FontSize',12,'FontWeight','bold')
  
  axes(ha(4))
  % subplot(224)
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,
    plot([L.X-sc*sind(L.FMSHmax),L.X+sc*sind(L.FMSHmax)]',[L.Y-sc*cosd(L.FMSHmax),L.Y+sc*cosd(L.FMSHmax)]','k','linewidth',2)
    scatter(L.X,L.Y,50,L.FMSHmax,'^','filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'d) Y&H13 SHmax at seismic stations','FontSize',12,'FontWeight','bold')

  figure(7),clf
  plot(B.SHmax,B.FMSHmax,'o',L.FastDirection,L.FMSHmax,'^','linewidth',1)
  axis equal,xlim([-90,90]),ylim([-90,90]),grid
  xticks((-90:30:90)),yticks((-90:30:90))
  legend('boreholes with breakouts','seismic stations with SWS fast directions')
  xlabel('observed SHmax or Fast Direction azimuth (degEofN)')
  ylabel('SHmax from Y&H13 regional focal mechanism inversion (degEofN)')
  title('clearly need to do this analysis considering local focal mechanisms...')

%
% Make an omni figure for the SCEC proposal, with all the relevant bits in the right arrangement
%

  figure(8),clf

  % set(gcf,'Position',[51,266,1103,690])
  % [ha,pos]=tight_subplot(3,2,0.05,0.15,0.05); % this just helps it be a nicer looking journal figure
 
  % axes(ha(1))
  subplot(321)
  % subplot(3,3,[1,2])
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,
    plot([B.X-sc*sind(B.SHmax),B.X+sc*sind(B.SHmax)]',[B.Y-sc*cosd(B.SHmax),B.Y+sc*cosd(B.SHmax)]','k','linewidth',2)
    scatter(B.X,B.Y,50,B.SHmax,'filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'a) borehole SHmax','FontSize',12,'FontWeight','bold')

  % axes(ha(3))
  subplot(323)
  % subplot(3,3,[4,5])
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,

    scatter(T.X,T.Y,5,T.fast,'filled'),

    plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',2)
    scatter(L.X,L.Y,50,L.FastDirection,'^','filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'b) SWS fast direction','FontSize',12,'FontWeight','bold')
  
  % axes(ha(5))
  subplot(325)
  % subplot(3,3,[7,8])
    plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok',L.X,L.Y,'^k'),
    axis equal,axis(P.R2)
    set(gca,'Color',[1 1 1]*0.75)
    hold on,
    plot([B.X-sc*sind(B.FMSHmax),B.X+sc*sind(B.FMSHmax)]',[B.Y-sc*cosd(B.FMSHmax),B.Y+sc*cosd(B.FMSHmax)]','k','linewidth',2)
    scatter(B.X,B.Y,50,B.FMSHmax,'filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    plot([L.X-sc*sind(L.FMSHmax),L.X+sc*sind(L.FMSHmax)]',[L.Y-sc*cosd(L.FMSHmax),L.Y+sc*cosd(L.FMSHmax)]','k','linewidth',2)
    scatter(L.X,L.Y,50,L.FMSHmax,'^','filled','markeredgecolor','k'),csym(90),colormap(cpolar),
    colorbar('YTick',(-90:30:90))
    xticks([]),yticks([])
    text(227,3737,'c) Y&H13 SHmax','FontSize',12,'FontWeight','bold')
  
  % [ha,pos]=tight_subplot(1,2,0.05,0.15,0.05); % this just helps it be a nicer looking journal figure

  % axes(ha(2))
  subplot(122)
  % subplot(133)
    plot(B.SHmax,B.FMSHmax,'o',L.FastDirection,L.FMSHmax,'^','linewidth',1)
    axis equal,xlim([-90,90]),ylim([-90,90]),grid
    xticks((-90:30:90)),yticks((-90:30:90))
    legend('boreholes with breakouts','seismic stations with SWS fast directions',...
      'location','northwest','FontSize',12,'FontWeight','bold')
    xlabel('observed BH SHmax or SWS Fast Direction (degEofN)','FontSize',12)
    ylabel('SHmax from Y&H13 regional focal mechanism inversion (degEofN)','FontSize',12)
    text(-88,-85,'d) ','FontSize',12,'FontWeight','bold')
    set(gcf,'renderer','Painters')

toc