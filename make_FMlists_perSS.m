clear,tic
addpath kluttrell_local 
setenv('PATH', [getenv('PATH') ':/usr/local/bin']); 

%
% February 2020: one more update, deep-only subsets
%   - add two new depth criteria: z > 3km and z> 5 km, labeled as "zmin"
%   - that way we're comparing independent and complimentary SHmax estimates
%
% January 2020: "final" update for paper
%   - using the new Table 1 list of Boreholes (excludes some and adds some 
%     relative to what we had before). Also fixes some prior location issues
%   - using extended YHS catalogs from 1981 - 2018, increasing catalog by 20%
%   - updated set of distance criteria we're interested in, with sensible upper and lower limits
%
% January 2019: update to calculate max distance for the 6 "closest" categories
%   - had previously done this ad hoc in the comments of plot_hyperlocalSHmaxCompare_v3.m
%   - codify it here, and save the results to a mat file for use in plot_hyperlocalSHmaxPDFs
%
% September 2018: update post SCEC-presented analysis
%   - add two new distance categories to run FM inversions over:
%      dmax = 35 km, and closest 20
%   - Jeanne will rerun all inversions and store the full PDF of the SHmax calculation
%   - Meanwhile, use this as an opportunity to clean up the script to a "polished" final form
%
% August 2018:
%   - For each BH SHmax observation within the LA basin...
%   - Identify "hyper local" subsets of nearby FMs for stress inversion
%      - Quality = Any (we can't be too picky, not enough events)
%      - Depth = Any, < 5 km (shallow), < 3 km (super shallow)... many of these will be empty
%      - dmax = closest 10 regardless of distance, < 1km, 2, 5, 10, 15, 20, 50
%      - time = Any
%   - Save the subsets to individual files that can be sent to Jeanne
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE DATA SETS
%
BHinfile='dats/Table_S1.csv';
FMinfile='dats/FM_subsets_asof_20200108.shmax';
SSinfile='dats/LiPeng2017tableS3.csv';

[B,S,L,F,Regions,P]=load_and_parse_SSBHFM(BHinfile,FMinfile,SSinfile);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the basic locations and subregions, just so we can tell what's up
%  - borehole depths are all < 5 km and most < 3 km, so these really are good 
%    categories to use for "shallow" vs "deep"
%
  sc=12;
  P.R4=[200 800 3500 4200];
  figure(1),clf
  plot(P.cutm(:,1),P.cutm(:,2),'b',P.futm(:,1),P.futm(:,2),'k'),axis equal,axis(P.R4)
  hold on,
  plot(F.X,F.Y,'.','color',[1 1 1]*0.75) % add the seimicity in the background
  plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',1)
  scatter(L.X,L.Y,30,L.ResultantLength,'filled'),colorbar,grid %csym(90)
  %yticks(3730:20:3810)
  %colormap(cpolar),
  %cbh=colorbar;set(cbh,'YTick',[-90:30:90])
  title('FastDirections for all the SS colored by ResultantLength') 
  

  %for k=1:length(Regions.R)
   % [Regions.boxX{k,1},Regions.boxY{k,1}]=ll2utm(Regions.boxy{k},Regions.boxx{k}); % utm for determining distances
   % Regions.boxX{k}=Regions.boxX{k}/1e3;Regions.boxY{k}=Regions.boxY{k}/1e3; % km are easier to deal with than m
   % plot(Regions.boxX{k},Regions.boxY{k},'k','linewidth',1)
  %end

  %figure(2),clf,
  %histogram(L.Nmeasurements),hold on,histogram(L.ResultantLength)
  %legend('deepest breakout','shallowest breakout'),xlabel('m'),ylabel('count')

  % pulling out just the LAbasin region
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
  
  figure(2),clf
  plot(P.cutm(:,1),P.cutm(:,2),'b',P.futm(:,1),P.futm(:,2),'k'),axis equal,axis(P.R4)
  hold on,
  plot(F.X,F.Y,'.','color',[1 1 1]*0.75) % add the seimicity in the background
  plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',1)
  scatter(L.X,L.Y,30,L.ResultantLength,'filled'),colorbar,grid %csym(90)
  title ('Only in LABasin')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT THROUGH AND IDENTIFY SUBSETS
%  - 57 BH x 3 z criteria x 14 distance criteria = 2394
%  - plus closest 20 and closest 10 = 2736
%
  zmaxset = [100 5 3]; % all, shallow, shallowest
  dmaxset = [20,15,10,8,6,4,2,1]; % far to near
%dmaxset = [35,30,25,20,18,16,14,12,10,8,6,4,2,1]; % far to near older
  
  toc
  for k=1:numel(L.X)
    D{k,1}=sqrt((F.X-L.X(k)).^2+(F.Y-L.Y(k)).^2); % distance from each FM to this BH

    for j=1:numel(dmaxset)
      Iall{k,j}=find(D{k}<dmaxset(j));
      Ishallow{k,j}=find(D{k}<dmaxset(j) & F.z < zmaxset(2));
      Ishallowest{k,j}=find(D{k}<dmaxset(j) & F.z < zmaxset(3));
      Ideep{k,j}=find(D{k}<dmaxset(j) & F.z >= zmaxset(2));
      Ideepest{k,j}=find(D{k}<dmaxset(j) & F.z >= zmaxset(3));
    end

    % sort the FMs in each depth class in preparation for finding the closest ones
    %[Dsorted,isorted]=sort(D{k});
    %isortedshallow=find(F.z(isorted)<zmaxset(2));
    %isortedshallowest=find(F.z(isorted)<zmaxset(3));
    %isorteddeep=find(F.z(isorted)>=zmaxset(2));
    %isorteddeepest=find(F.z(isorted)>=zmaxset(3));

    % find the 10 closest ones in each depth class
    %j=j+1;
    %Nclose=10;
    %Iall{k,j}=isorted(1:Nclose);
    %Ishallow{k,j}=isorted(isortedshallow(1:Nclose));
    %Ishallowest{k,j}=isorted(isortedshallowest(1:Nclose));
    %Ideep{k,j}=isorted(isorteddeep(1:Nclose));
    %Ideepest{k,j}=isorted(isorteddeepest(1:Nclose));

    % find the 20 closest ones in each depth class
    %j=j+1;
    %Nclose=20;
    %Iall{k,j}=isorted(1:Nclose);
    %Ishallow{k,j}=isorted(isortedshallow(1:Nclose));
    %Ishallowest{k,j}=isorted(isortedshallowest(1:Nclose));
    %Ideep{k,j}=isorted(isorteddeep(1:Nclose));
    %Ideepest{k,j}=isorted(isorteddeepest(1:Nclose));
  end
  toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW MANY OF THESE ARE EMPTY?  HOW FAR OUT IS THE 10th/20th CLOSEST EQ?
%
  Nall = cell2mat(cellfun(@numel,Iall,'uniformoutput',false));
  Nshallow = cell2mat(cellfun(@numel,Ishallow,'uniformoutput',false));
  Nshallowest = cell2mat(cellfun(@numel,Ishallowest,'uniformoutput',false));
  Ndeep = cell2mat(cellfun(@numel,Ideep,'uniformoutput',false));
  Ndeepest = cell2mat(cellfun(@numel,Ideepest,'uniformoutput',false));

  Nempty=numel(find(Nall==0))+numel(find(Nshallow==0))+numel(find(Nshallowest==0))+numel(find(Ndeep==0))+numel(find(Ndeepest==0));
  Nsmall=numel(find(Nall<6))+numel(find(Nshallow<6))+numel(find(Nshallowest<6))+numel(find(Ndeep<6))+numel(find(Ndeepest<6));

  %h=@(x,y)x(y); % this is the function we're going to use to pull the correct value out of the distance from BH cell array

  %idmaxofeach=cellfun(@max,Iall(:,end-1),'uniformoutput',false); %end-1 gets you to the "10 closest" column
  %d10th=cellfun(h,D,idmaxofeach); % distance to 10th closest FM
  %idmaxofeach=cellfun(@max,Ishallow(:,end-1),'uniformoutput',false); % same for z<5km
  %d10th_shallow=cellfun(h,D,idmaxofeach);
  %idmaxofeach=cellfun(@max,Ishallowest(:,end-1),'uniformoutput',false); % same for z<3km
  %d10th_shallowest=cellfun(h,D,idmaxofeach);                         
%   idmaxofeach=cellfun(@max,Ideep(:,end-1),'uniformoutput',false); % same for z>=5km
%   d10th_deep=cellfun(h,D,idmaxofeach);
%   idmaxofeach=cellfun(@max,Ideepest(:,end-1),'uniformoutput',false); % same for z>=3km
%   d10th_deepest=cellfun(h,D,idmaxofeach);                         
% 
%   idmaxofeach=cellfun(@max,Iall(:,end),'uniformoutput',false);
%   d20th=cellfun(h,D,idmaxofeach); % distance to 20th closest FM
%   idmaxofeach=cellfun(@max,Ishallow(:,end),'uniformoutput',false); % same for z<5km
%   d20th_shallow=cellfun(h,D,idmaxofeach);
%   idmaxofeach=cellfun(@max,Ishallowest(:,end),'uniformoutput',false); % same for z<3km
%   d20th_shallowest=cellfun(h,D,idmaxofeach);                         
%   idmaxofeach=cellfun(@max,Ideep(:,end),'uniformoutput',false); % same for z>=5km
%   d20th_deep=cellfun(h,D,idmaxofeach);
%   idmaxofeach=cellfun(@max,Ideepest(:,end),'uniformoutput',false); % same for z>=3km
%   d20th_deepest=cellfun(h,D,idmaxofeach);                         

  % This file is read in by plot_hyperlocalSHmaxPDFs.m, in order to color FM SHmax fits by distance to 10th/20th EQ 
  % save('dats/hyperlocal_NEQ_and_dmax_perCriteria','Nall','Nshallow','Nshallowest','d10th','d10th_shallow','d10th_shallowest','d20th','d20th_shallow','d20th_shallowest')


  figure(3),clf,CA=[-0.5,4];
  subplot(151),imagesc(log10(Nall)),colorbar,caxis(CA),title('log10(Nall)'),ylabel('Seismic Station #')
    xticks(1:numel(dmaxset)+2),xticklabels([string(dmaxset)]),xlabel('dmax (km)')
  subplot(152),imagesc(log10(Nshallow)),colorbar,caxis(CA),title('log10(Nshallow)')
    xticks(1:numel(dmaxset)+2),xticklabels([string(dmaxset)]),xlabel('dmax (km)')
  subplot(153),imagesc(log10(Nshallowest)),colorbar,caxis(CA),title('log10(Nshallowest)')
    xticks(1:numel(dmaxset)+2),xticklabels([string(dmaxset)]),xlabel('dmax (km)')
  subplot(154),imagesc(log10(Ndeep)),colorbar,caxis(CA),title('log10(Ndeep)')
    xticks(1:numel(dmaxset)+2),xticklabels([string(dmaxset)]),xlabel('dmax (km)')
  subplot(155),imagesc(log10(Ndeepest)),colorbar,caxis(CA),title('log10(Ndeepest)')
    xticks(1:numel(dmaxset)+2),xticklabels([string(dmaxset)]),xlabel('dmax (km)')
  %subplot(166),plot(d10th,1:numel(L.X),'.-',d20th,1:numel(L.X),'.-'),axis ij,title('distance to closest'),xlabel('distance (km)'),grid,legend('10th','20th','location','southeast')

toc

%stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that I have my lists of mechanisms around each borehole,
%  - use the indices to sed out these lines from th YSH catalog and print them out to individual files.
%    note: this is inefficient... will takes an hour to do all of them this way...
%  - alternatively, write out the lines using matlab (already read them in)
%    this is ~100x faster
%  - takes a few minutes
%
  %YSH_without_Q=cell2mat(C(1:20));

  % DIR='dats/FM_subsets_asof_20200108/'; % these were the all, shallow, and shallowest ones
  DIR='dats/FM_subsets_asof_20210615/'; % these are the seismic station lists

  tic
  for k=1:numel(L.X)
    %
    % do the first set of distance criteria, the ones with a defined dmax.
    %
      for j=1:numel(dmaxset)
        filename=[DIR,'SS_',L.station{k},'_dmax',num2str(dmaxset(j),'%02d'),'km_zmaxany.dat'];
          fid=fopen('temp1','w');
          fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Iall{k,j},:)');
          fclose(fid);

          fid=fopen(['temp2'],'w');
          fprintf(fid,'%s\n',string(F.Q(Iall{k,j})));
          fclose(fid);

          system(['paste temp1 temp2 > ',filename]);
          system('rm temp1 temp2');

        filename=[DIR,'SS_',L.station{k},'_dmax',num2str(dmaxset(j),'%02d'),'km_zmax5km.dat'];
          fid=fopen('temp1','w');
          fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ishallow{k,j},:)');
          fclose(fid);

          fid=fopen(['temp2'],'w');
          fprintf(fid,'%s\n',string(F.Q(Ishallow{k,j})));
          fclose(fid);

          system(['paste temp1 temp2 > ',filename]);
          system('rm temp1 temp2');

        filename=[DIR,'SS_',L.station{k},'_dmax',num2str(dmaxset(j),'%02d'),'km_zmax3km.dat'];
          fid=fopen('temp1','w');
          fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ishallowest{k,j},:)');
          fclose(fid);

          fid=fopen(['temp2'],'w');
          fprintf(fid,'%s\n',string(F.Q(Ishallowest{k,j})));
          fclose(fid);

          system(['paste temp1 temp2 > ',filename]);
          system('rm temp1 temp2');

        filename=[DIR,'SS_',L.station{k},'_dmax',num2str(dmaxset(j),'%02d'),'km_zmin5km.dat'];
          fid=fopen('temp1','w');
          fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ideep{k,j},:)');
          fclose(fid);

          fid=fopen(['temp2'],'w');
          fprintf(fid,'%s\n',string(F.Q(Ideep{k,j})));
          fclose(fid);

          system(['paste temp1 temp2 > ',filename]);
          system('rm temp1 temp2');

        filename=[DIR,'SS_',L.station{k},'_dmax',num2str(dmaxset(j),'%02d'),'km_zmin3km.dat'];
          fid=fopen('temp1','w');
          fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ideepest{k,j},:)');
          fclose(fid);

          fid=fopen(['temp2'],'w');
          fprintf(fid,'%s\n',string(F.Q(Ideepest{k,j})));
          fclose(fid);

          system(['paste temp1 temp2 > ',filename]);
          system('rm temp1 temp2');
      end

    %
    % do the n-1th distance constraint: the closest 10 has a different naming convention
    %
%       j=j+1; 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest10_zmaxany.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Iall{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Iall{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest10_zmax5km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ishallow{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ishallow{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest10_zmax3km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ishallowest{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ishallowest{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest10_zmin5km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ideep{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ideep{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest10_zmin3km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ideepest{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ideepest{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%     %
%     % do the nth distance constraint: the closest 20 has a different naming convention
%     %
%       j=j+1;
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest20_zmaxany.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Iall{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Iall{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest20_zmax5km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ishallow{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ishallow{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest20_zmax3km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ishallowest{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ishallowest{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest20_zmin5km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ideep{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ideep{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');
% 
%       filename=[DIR,'SS_',L.station{k},num2str(dmaxset(j),'%02d'),'_dmaxClosest20_zmin3km.dat'];
%         fid=fopen('temp1','w');
%         fprintf(fid,'%4d %2d %2d %2d %2d %6.3f %8d %9.5f %10.5f %7.3f %6.3f %4d %3d %4d %3d %3d %4d %5.2f %4d %5.2f\n',F.YSH_without_Q(Ideepest{k,j},:)');
%         fclose(fid);
% 
%         fid=fopen(['temp2'],'w');
%         fprintf(fid,'%s\n',string(C{21}(Ideepest{k,j})));
%         fclose(fid);
% 
%         system(['paste temp1 temp2 > ',filename]);
%         system('rm temp1 temp2');

    toc
   
  end
stop



    

