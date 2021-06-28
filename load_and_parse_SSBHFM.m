function [B,S,L,F,Regions,P]=load_and_parse_SSBHFM(BHinfile,FMinfile,SSinfile)
%
% VERY IMPORTANT COMMENT, REGARDLESS OF SPELLING!!!!!!!
% EXTRAORDINARILY IMPORTANT!
%
%NEW NEW
%
% [B,S,F,Regions,P]=load_and_parse_BHFM(BHinfile,FMinfile)
%
% inputs: 
%  - BHinfile to load and parse, e.g. dats/Table_S1.csv
%  - FMinfile to load and parse, e.g. dats/dats/FM_subsets_asof_20200108.shmax
%  - SSinfile to load and parse, e.g. dats/LiPeng2017tableS3.csv
% outputs:
%  - B structural array with BH info
%  - S structural array with FM inversion info
%  - L structural array with SWS info 
%  - F structural array with inividual focal mechanism info, from dats/YHS_fm_1981_2018.hash
%  - Regions structural array with geometry info, from dats/subregionlist.csv
%  - P structural array with plotting info from 
%      - dats/SoCal_coastfile.xy
%      - dats/allCAfaults.MATLAB.dat
%      - grds/Yang.FMSHmax.grd
%      - grds/CVM15_10sec.grd
%      - regional bounds for geographic and utm coordinates
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME THE FILES WE'LL BE LOADING IN BELOW (so they're easily changed later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % clear,tic
  BHinfile='dats/Table_S1.csv';
  FMinfile='dats/FM_subsets_asof_20200108.shmax';
  SSinfile='dats/LiPeng2017tableS3.csv';

  ZSTATSfile='dats/zstats_for_FM_subsets_asof_20200108.txt';
  NEQfile='dats/FM_subsets_asof_20200108.Neqdmax.mat';
  HASHfile='dats/YHS_fm_1981_2018.hash';
  REGIONSfile='dats/subregionlist.csv';
  COASTfile='dats/SoCal_coastfile.xy';
  FAULTfile='dats/allCAfaults.MATLAB.dat';
  YH13grdfile='grds/Yang.FMSHmax.grd';
  CVMgrdfile='grds/CVM15_10sec.grd';
  CVMgrdfileUTM='grds/CVM15_250m_UTM11.grd';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE BOREHOLE DATA AND REGIONAL FM INVERSION DATA
% This represents the state of knowledge before we came along...
%  - borehole breakout info: from January 2020 update
%  - regional inversions from Yang 2013
%  - Also comparison info: Community Velocity Model depth to basement, etc.
%  - and map info for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % 1) Load the Seismic Stations 
  %
  
      fid=fopen(SSinfile);
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
     
    L.FastDirection(find(L.FastDirection > 90))=L.FastDirection(find(L.FastDirection > 90))-180; % make boreholes [-90,90] to match FM SHmax
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % 1b) Load the borehole SHmax observations
  %
    % fid=fopen('dats/Table_S1.csv');
    fid=fopen(BHinfile);
    C=textscan(fid,'%f %f %f %f %f %f %f %s %f %s %s %s %s %s %f %s %f %f %s %f %f','HeaderLines',1,'delimiter',',');
    fclose(fid);

    B.x=C{2};
    B.y=C{1};
    B.z1=C{3}; % top of breakouts
    B.z2=C{4}; % bottom of breakouts
    B.SHmax=C{5}; % azimuth of breakouts (degEofN)
    B.stdSHmax=C{6}; % std of azimuth of breakouts (deg)
    B.dSHmax=C{7}; % functional uncertainty of azimuth of breakouts (deg) (larger of std or 10º)
    % B.dSHmax=max(C{6},15); % try larger of std or 15º

    B.FMSHmax=C{9}; % SHmax from Y&H2013
    B.Ref=C{10}; % which paper is this from?
    B.REG=C{16}; % Sub Region

    [B.X,B.Y]=ll2utm(B.y,B.x); % utm for determining distances
    B.X=B.X/1e3;B.Y=B.Y/1e3; % km are easier to deal with than m

    B.SHmax(find(B.SHmax > 90))=B.SHmax(find(B.SHmax > 90))-180; % make boreholes [-90,90] to match FM SHmax
    
    %
    % FM SHmax UNCERTAINTIES from regionally smoothed model:
    %  - Y&H2013 report S1 azimuth uncertainty of 11º across most of the modeled region
    %
      B.FMdSHmax=ones(size(B.FMSHmax))*11;

    %
    % Define reference number for plotting
    %
      [~,~,B.RefNum]=unique(B.Ref);

    %
    % Order number in Table_S1_v1.csv.  This was the order used for the FM subsets originally,
    % so it essentially serves as a unique id# for each BH.  But we want to sort the BHs more
    % sensibly based on our results for display in the journal table, so will need to reorder
    % the S.* results below to follow the new order.  Do this within the main program, to
    % preserve backwards compatibility with the older stacklocalPDFs_v1.m
    %
      B.id=C{18};

    %
    % Geologic Context, string and number
    %
      B.Geol=C{19};
      B.GeolNum=C{20};
      % [~,~,B.GeolNum]=unique(B.Geol,'stable');
      
  % 
  % 2) Load the FM catalog, and the Y&H inversion full grid 
  %    (values sampled at BH locations loaded above)
  %
    % fid=fopen('dats/YHS_fm_1981_2018.hash');
    fid=fopen(HASHfile);
    C=textscan(fid,[repmat('%f ',1,20),'%s']);
    fclose(fid);

    F.t=datenum(cell2mat(C(1:6)));
    F.x=C{9};
    F.y=C{8};
    F.z=C{10};
    F.M=C{11};
    F.Q=C{21}; %mechanism quality...
    [F.X,F.Y]=ll2utm(F.y,F.x); % utm for determining distances
    F.X=F.X/1e3;F.Y=F.Y/1e3; % km are easier to deal with than m
    F.YSH_without_Q=cell2mat(C(1:20)); 

    % 
    % Load the Yang & Hauksson 2013 inversion results, the full grid
    %
      % [P.YH13.x,P.YH13.y,P.YH13.z]=grdread('grds/Yang.FMSHmax.grd');
      [P.YH13.x,P.YH13.y,P.YH13.z]=grdread(YH13grdfile);

  % 
  % 3) Load the Community Velocity Model depth to basement
  %  - do the 10sec gridded version of the top surface,
  %    which will be high enough resolution for our purposes
  %  - depth to shallowest basement surface, in negative meters
  %
    % [P.CVM.x,P.CVM.y,P.CVM.z]=grdread('grds/CVM15_10sec.grd');
    [P.CVM.x,P.CVM.y,P.CVM.z]=grdread(CVMgrdfile);
    [P.CVM.X,P.CVM.Y,P.CVM.Z]=grdread(CVMgrdfileUTM);
   
    [~,r]=system('cat dats/Table_S1.csv | sed ''1d'' | awk -F, ''{print $2,$1}'' | gmt grdtrack -Ggrds/CVM15_10sec.grd -Z');
    C=textscan(r,'%f');
    L.CVMz=C{1};

  %
  % 3.5) Determine the abosolute elevations of breakouts
  %  - track elevation at each point from topo grd
  %  - subtract depths from surface elevation to get breakout elevations
  %
    [~,r]=system('cat dats/Table_S1.csv | sed ''1d'' | awk -F, ''{print $2,$1}'' | gmt grdtrack -Ggrds/topo.grd -Z');
    C=textscan(r,'%f');
    B.topoz=C{1};
    B.z1el=B.topoz-B.z1;
    B.z2el=B.topoz-B.z2;


  %
  % 4) Load info for plotting
  %  - faults, coastline, regional bounds
  %  - define subregions (these are based off the subsets inherited from W&S)
  %
    %
    % coastline, faults, and regional bounds for plotting
    %
      % P.c=load('dats/SoCal_coastfile.xy');
      % P.f=load('dats/allCAfaults.MATLAB.dat');
      P.c=load(COASTfile);
      P.f=load(FAULTfile);
      P.cutm=ll2utm(P.c(:,2),P.c(:,1),11);P.cutm=P.cutm/1e3;
      P.futm=ll2utm(P.f(:,2),P.f(:,1),11);P.futm=P.futm/1e3;

      P.R=[-120,-117.8,33.7,34.5]; % LA basin focus
      P.R2=[225,430,3730,3825]; % LA basin focus in UTM km
      P.R3=[-120,-117.76,33.7,34.54]; % best match to UTM

    %
    % Define 8 subregions based on the locations of the borehole observations
    %  - pull it from the defined file, but that file was made from the subregions
    %    defined in an earlier version of this list
    %
      %
      % load the basic subregion info
      %
        % fid=fopen('dats/subregionlist.csv');
        fid=fopen(REGIONSfile);
        C=textscan(fid,'%f %f %f %f %s','delimiter',',');
        fclose(fid);

        Regions.name=C{5};
        Regions.R=mat2cell(cell2mat(C(1:4)),ones(1,numel(C{1})));

    %
    % define the vectors that would draw a box around each subregion, both in lat/lon and in utm km
    %
      for k=1:length(Regions.R)
        Regions.boxx{k,1}=Regions.R{k}([1 2 2 1 1])';
        Regions.boxy{k,1}=Regions.R{k}([3 3 4 4 3])';
        [Regions.boxX{k,1},Regions.boxY{k,1}]=ll2utm(Regions.boxy{k},Regions.boxx{k}); % utm for determining distances
        Regions.boxX{k}=Regions.boxX{k}/1e3;Regions.boxY{k}=Regions.boxY{k}/1e3; % km are easier to deal with than m
      end

    %
    % Which borehole points are in which subregions?
    %  - there's now a column about this in Table 1, which was derived from a previous version of this script
    %
      for k=1:length(Regions.R)
        L.RegionNum(find(strcmp(L.Nmeasurements,Regions.name{k})),1)=k;
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD HYPERLOCAL FM INVERSIONS FOR SHmax
% THESE ARE THE RESULTS AS I GET THEM FROM JEANNE
%  - this time, Jeanne gave me her best SHmax and also the PDF of SHmax uncertainty
%  - Recall, this is because all we have from BH is SHmax, so that's what we have to compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % %
  % % Load the results from localized inversions
  % %  - instead of S123 trend/plunge and uncertainty,
  % %    this time we have Jeanne's estimate of best SHmax, and the full PDF of SHmax values
  % %  - So, I shouldn't have to fret about the best way to deal with SHmax uncertainties or be sure I'm doing it right
  % %  - NOTE: this long list doesn't really get used in the calling scripts (only the parsed grids), so load these
  % %    variables into a different structure than the one that will ultimately be passed.
  % %
  %   FILENAME='dats/FM_subsets_asof_20200108.shmax';
     fid=fopen(FMinfile);
     C=textscan(fid,['%f %s',repmat(' %f',1,38)]);
     fclose(fid);

     R.N1=C{1}; % file number, in alphabetical order?
     R.filename=C{2};
     R.N2=C{3}; % seems to be the same as N1, PROBABLY AN INDEX FOR PASTING THESE TWO SECTIONS TOGETHER
     R.SHmax=C{4};
     R.SHmax_pdf=cell2mat(C(5:40));

  % %
  % % parse the info in the filename to figure out which borehole and criteria each result belongs to
  % %
  %   %
  %   % borehole number, 1-60
  %   %
       [~,R.BH]=system(['awk ''{print $2}'' ',FMinfile,' |  awk -F/ ''{print $2}'' | awk -F_ ''{print $1}'' | awk -FH ''{print $2}'' ']);
       R.BH=str2num(R.BH);

  %   %
  %   % max depth: all (999 km), shallow (5 km), or super shallow (3 km)
  %   %
       [~,R.zmax]=system(['awk ''{print $2}'' ',FMinfile,' | awk -F/ ''{print $2}'' | awk -F_ ''{print $3}'' | sed ''s/any/999km/'' | awk -Fx ''{print $2}'' | awk -Fk ''{print $1}''']);
       R.zmax=str2num(R.zmax);

  %   %
  %   % max distance: # of km, -10 for closest 10 events, -20 for closest 20 events
  %   %
       [~,R.dmax]=system(['awk ''{print $2}'' ',FMinfile,' | awk -F/ ''{print $2}'' | awk -F_ ''{print $2}'' | sed ''s/Closest10/-10km/'' | sed ''s/Closest20/-20km/'' | awk -Fx ''{print $2}'' | awk -Fk ''{print $1}''']);
       R.dmax=str2num(R.dmax);

  %   %
  %   % Based on Borehole Number, assign the appropriate locations...
  %   %
       R.x=B.x(R.BH);
       R.y=B.y(R.BH);
       R.X=B.X(R.BH);
       R.Y=B.Y(R.BH);
       R.z1=B.z1(R.BH);
       R.z2=B.z2(R.BH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD depth statistics per FM solution (calculated by calcFMdepthstats.com)
%  - to be parsed into 3-D arrays below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fid=fopen(ZSTATSfile);
   C=textscan(fid,'BH%f_dmax%fkm_zmax%fkm.dat %f %f %f %f %f %f');
   fclose(fid);

   rBH=C{1};       % borehole number
   rdmax=C{2};     % dmax in km (negative for "closest")
   rzmax=C{3};     % zmax in km (999 km for "any")
   rzmean=C{4};    % mean depth of eqs in this set, NaN if empty
   rz05th=C{5};    % 5th depth quantile (so, top) of eqs in this set, NaN if empty
   rz25th=C{6};    % 25th depth quantile of eqs in this set, NaN if empty
   rz50th=C{7};    % median depth of eqs in this set, NaN if empty
   rz75th=C{8};    % 75th depth quantile of eqs in this set, NaN if empty
   rz95th=C{9};    % 95th depth quantile (so, base) of eqs in this set, NaN if empty

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARSE OUT HYPERLOCAL SHmax VALUES INTO SENSIBLE 3D GRID
%  - for each BH, zmax, and dmax option.  Many of these will be empty.
%  - 3D array with dimensions of NBH x Ndmax x Nzmax
%  - While we're at it, use the same loops to parse the FM depth statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %
    % % Define parameter sets used
    % %   - list of zmax, dmax, and BH# considered
    % %   - Bins associted with the SHmax pdf distribution
    % %
       S.zmaxset=unique(R.zmax);
       S.dmaxset=unique(R.dmax);
       S.BHset=[1:numel(B.x)]';
       S.bins=-87.5:5:87.5;
        
    % %
    % % Define the empty 3-D arrays to be filled
    % %

       S.SHmax=ones(numel(S.BHset),numel(S.dmaxset),numel(S.zmaxset))*NaN;

    % %
    % % loop to fill in the arrays with available data
    % %
       for i=1:numel(S.BHset)
        for j=1:numel(S.dmaxset)
          for k=1:numel(S.zmaxset)

            id=find(R.BH==S.BHset(i) & R.dmax==S.dmaxset(j) & R.zmax==S.zmaxset(k));
            if numel(id)>0
              S.SHmax(i,j,k)=R.SHmax(id); % best SHmax, Jeanne's estimate
              S.SHmax_pdf{i,j,k}=R.SHmax_pdf(id,:); % pdf of SHmax, Jeanne's estimate
            end

            id2=find(rBH==S.BHset(i) & rdmax==S.dmaxset(j) & rzmax==S.zmaxset(k));

            S.zmean(i,j,k)=rzmean(id2);
            S.z05th(i,j,k)=rz05th(id2);
            S.z25th(i,j,k)=rz25th(id2);
            S.z50th(i,j,k)=rz50th(id2);
            S.z75th(i,j,k)=rz75th(id2);
            S.z95th(i,j,k)=rz95th(id2);
          end
        end
      end

    % %
    % % Make comparably sized arrays of Borehole SHmax and dSHmax, as well as regional FM SHmax and dSHmax
    % %
       S.B_SHmax=repmat(B.SHmax,[1,numel(S.dmaxset),numel(S.zmaxset)]);
       S.B_dSHmax=repmat(B.dSHmax,[1,numel(S.dmaxset),numel(S.zmaxset)]);
       S.F_SHmax=repmat(B.FMSHmax,[1,numel(S.dmaxset),numel(S.zmaxset)]);
       S.F_dSHmax=repmat(B.FMdSHmax,[1,numel(S.dmaxset),numel(S.zmaxset)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD Number of earthquakes per solution and distance to 10th or 20th closest earthquake
%  - parse out number of earthquakes per set (note: the dmaxset list is in the opposite order from what
%    is established by reading in the Jeanne results above, hence the need for the fliplr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   T=load(NEQfile);
   S.Neq(:,:,1)=fliplr(T.Nshallowest);
   S.Neq(:,:,2)=fliplr(T.Nshallow);
   S.Neq(:,:,3)=fliplr(T.Nall);
   S.d10th=[T.d10th_shallowest,T.d10th_shallow,T.d10th];
   S.d20th=[T.d20th_shallowest,T.d20th_shallow,T.d20th];

  %S=[]; % dummy assign to keep it from complaining...
% toc
