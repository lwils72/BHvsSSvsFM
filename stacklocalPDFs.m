clear
addpath kluttrell_local 
setenv('PATH', [getenv('PATH') ':/usr/local/bin']);

%
% stacklocalPDF.m: See previous stacklocalPDF, version 1
%  - Analysis was completed in previous version 1.
%  - This version aims to tighted figures for journal presentation, rather than complete analysis.
%  - RECALL: in version 1, we made all analaysis figures both using MAD and RMS metrics of stacking
%    - for MAD, we did FFT convolution to "analytically" calculate the full pdf
%    - for RMS, draw lots of samples from the difference PDF (still calculated "analytically" using
%      FFTs and the convolution theorem), then square, sum, divide by N, and sqrt: more of a bootstrap
%  - For this version, we're dropping the MAD ones and going instead with the RMS ones
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD AND PLOT THE STUFF WE ALREADY KNOW ABOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Load the data and info we need, and define functions for comparison metrics
  %
    tic
    BHinfile='dats/Table_S1.csv';
    FMinfile='dats/FM_subsets_asof_20210615.shmax';
    SSinfile='dats/LiPeng2017tableS3.csv';
    [B,S,L,F,Regions,P]=load_and_parse_SSBHFM(BHinfile,FMinfile,SSinfile); %make sure everything is in correct order with infiles uncommented

   %
   % This excludes SS that are not in the LA basin 
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

    % THESE ACTUALLY AREN"T BEING USED THAT MUCH... CAN... UM... clear these a bit?
    daz = @(x,y) (mod(y-x+90,180)-90);  % Circular Difference
    rmsaz = @(x,y) sqrt(sum(daz(x,y).^2,'omitnan')./sum(~isnan(y))); % RMS of Circular Difference, ignores NaNs, works with arrays of column vectors
    toc
  %
  % Mapview plot of borehole locations, colored by number because it helps to easily identify them.
  %
    figure(1),clf
    subplot(211)
    plot(P.cutm(:,1),P.cutm(:,2),'k'),hold on,
    plot(P.futm(:,1),P.futm(:,2),'color',[1 1 1]*0.75)
    axis equal,axis(P.R2)
    scatter3(B.X,B.Y,1:numel(B.X),50,1:numel(B.X),'filled'),colorbar
    plot3(B.X,B.Y,1:numel(B.X),'ok')
    colormap(jet),view(2) % sets view to X-Y 2D 
    set(gcf,'renderer','Painters')

  %
  % depth distributions of FMs, as well as CVM depth and Breakout depths
  %  - oh, drat, don't have these as distributions, just the stats... oh well
  %
      %
      % fine earthquakes within the region and their associated histogram
      %  - can we do this in real time?
      %
        kmbinedges=-1:0.1:30; % For consistent depth bins

        iin=find(F.X>P.R2(1) & F.X<P.R2(2) & F.Y>P.R2(3) & F.Y<P.R2(4));
        [Neqz,edges]=histcounts(F.z(iin),kmbinedges,'normalization','pdf');
        zbins=edges(1:end-1)+diff(edges)/2;

        Ncvmz=histcounts(-P.CVM.Z(:)/1e3,kmbinedges,'normalization','pdf');

        Bzset=cell2mat(cellfun(@(x,y)x:y,num2cell(floor(B.z2el)),num2cell(ceil(B.z1el)),'uniformoutput',false)')';
        Nbhz=histcounts(-Bzset/1e3,kmbinedges,'normalization','pdf');

      subplot(2,1,2)
        SCALE=4;
        plot(1+Neqz*SCALE,-zbins,'k',2+Ncvmz*SCALE,-zbins,'k',3+Nbhz*SCALE/2,-zbins,'k')
        xlim([0.5,4]),xticks(1:3),xticklabels({'EQs','CVM','breakouts'})
        ylim([-20,1]),ylabel('depth (km)'),grid
        title('depth distributions within this area')

        hold on
        polyy=-[zbins(1),zbins,zbins(end)];
        polyx=[0,Neqz,0]*SCALE+1;
        pgon=polyshape(polyx,polyy);
        plot(pgon,'linestyle','none','facecolor',[1 1 1]*0.75)
        polyx=[0,Ncvmz,0]*SCALE+2;
        pgon=polyshape(polyx,polyy);
        plot(pgon,'linestyle','none','facecolor',[1 1 1]*0.75)
        polyx=[0,Nbhz,0]*SCALE/2+3;
        pgon=polyshape(polyx,polyy);
        plot(pgon,'linestyle','none','facecolor',[1 1 1]*0.75)

  %
  % Plot the 3D view of FMs, EQs, and CVM to scale (see previous figure 18)
  %
  L.Z=L.X*0;
  
    figure(2),clf
    scatter3(F.X,F.Y,-F.z,15,-F.z,'filled'),colorbar
    hold on,
    plot3(L.X,L.Y,L.Z/1e3,'ko','linewidth',2)
    plot3((B.X*[1 1])',(B.Y*[1 1])',[B.z1el B.z2el]'/1e3,'-r','linewidth',3)
    surf(P.CVM.X(1:2:end)/1e3,P.CVM.Y(1:2:end)/1e3,P.CVM.Z(1:2:end,1:2:end)/1e3,'facecolor','none') 
    plot(P.cutm(:,1),P.cutm(:,2),'k')
    axis equal,axis(P.R2),zlim([-30,1]),caxis([-20,0])
    %plot3([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]',[L.Z-sc*sind(L.FastDirection),L.Z+sc*sind(L.FastDirection)]','k','linewidth',0.05) %SS
    scatter3(L.X,L.Y,L.Z,10,L.ResultantLength,'filled'),colorbar,grid %csym(90)

    iNorthRidgeEQ=find(floor(F.t)==datenum([1994 01 17]) & F.M>6);
    plot3(F.X(iNorthRidgeEQ),F.Y(iNorthRidgeEQ),-F.z(iNorthRidgeEQ),'k*','markersize',20,'linewidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SORT THE FM INVERSION SOLUTIONS INTO NEW ORDER TO MATCH THE CURRENT BH LIST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Basically, anything within S.* that has 57 rows, the rows need to be reordered based on
  % the BH.id order.
  %  - EXCEPT!!! S.B_SHmax, because that was derived from B, so it already has the right order
  %
%          S.SSset=S.SSset(L.id,:);            % [57x1 double]
%          S.SHmax=S.SHmax(L.id,:,:);          % [57x16x3 double]
%          S.SHmax_pdf=S.SHmax_pdf(L.id,:,:);  % {57x16x3 cell}
%          S.zmean=S.zmean(L.id,:,:);          % [57x16x3 double]
%          S.z05th=S.z05th(L.id,:,:);          % [57x16x3 double]
%          S.z25th=S.z25th(L.id,:,:);          % [57x16x3 double]
%          S.z50th=S.z50th(L.id,:,:);          % [57x16x3 double]
%          S.z75th=S.z75th(L.id,:,:);          % [57x16x3 double]
%          S.z95th=S.z95th(L.id,:,:);          % [57x16x3 double]
%          S.B_SHmax=S.B_SHmax(B.id,:,:);      % [57x16x3 double] out on l&p
%          S.B_dSHmax=S.B_dSHmax(B.id,:,:);     % [57x16x3 double]
%          S.L_SHmax=S.F_SHmax(L.id,:,:);        % [57x16x3 double]
%          S.L_dSHmax=S.F_dSHmax(L.id,:,:);       % [57x16x3 double]
%          S.Neq=S.Neq(B.id,:,:);            % [57x16x3 double]
%          S.d10th=S.d10th(B.id,:);            % [57x3 double]
%          S.d20th=S.d20th(B.id,:);            % [57x3 double]

  %
  % Do the comparison and plot the raw results: 
  %  - This is actually only in the stacked plot?.  
  %  - we're overall much more interested in the full pdfs with uncertainty
  %  - DO THIS AFTER RE-SORTING S.* to avoid getting the wrong answers
  %
    %zBHvFM=abs(daz(S.SHmax,S.B_SHmax)); % Absolute Circular Difference of the Data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE MAPVIEW OF FM SOLUTIONS, JUST BECAUSE... THEY GO IN THE SUPPLEMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Test the Cheesy attempt to get the topography to show up in the background? (I'm so lazy...)
  %  - the GMT commands make just the gray shaded topo, like in actual Figure 1,
  %    which I then screen-shotted to make it into a png... did I mention I was lazy?
  %
    % gmt grdgradient grds/topo.grd -Nt.6 -A300 -Gtempgrad.grd
    % gmt makecpt -Cgray -T-26.5/1/.1  -Z   >  gr.cpt 
    % gmt psbasemap -JM5 -R-120/-117.75/33.7/34.55 -Ba0.2f0.2g0/a0.2f0.2g0":.""topo only":WSen -K -X1 -Y6 > figs/topoonly.ps
    % gmt psbasemap -JM5 -R-120/-117.75/33.7/34.55 -Ba0.2f0.2g0/a0.2f0.2g0":.""topo only":WSen -K -X1 -Y6 > figs/topoonly.ps
    % gmt grdimage tempgrad.grd -Itempgrad.grd -Cgr.cpt -JM -R -O -K >> figs/topoonly.ps
    I = imread('figs/topoonly.png'); 

    RI=[-120,-117.75,33.7,34.55];
    shift=[0.8,-1.4,-1.3,2.0]; % offset fudge factors due to projection wonkiness

    boxI=ll2utm([RI(3),RI(1);RI(3),RI(2);RI(4),RI(2);RI(4),RI(1)])/1e3;
    RIutm=[mean([boxI(1,1),boxI(4,1)]),mean(boxI(2:3,1)),mean(boxI(1:2,2)),mean(boxI(3:4,2))]+shift;

    % figure(99),clf,
    % image(RIutm(1:2),RIutm(3:4),flipud(I)),axis xy,
    % axis equal,axis tight
    % hold on
    % plot(P.cutm(:,1),P.cutm(:,2),'k',B.X,B.Y,'ok')
    % axis(RIutm)
    % grid

  %
  % This is the addition of the scatter plot from test_SWSvsFM figure (8) 
  % Different subplots in one large plot Zmax 5,10,any vs Dmax etc. 
  % Figures need to show the various z or d max distances. 
  
      figure(8),clf
      for j=1:numel(S.dmaxset)
         for k=1:numel(S.zmaxset)
      subplot(3,8,(k-1)*8+j)
      plot(L.FastDirection,S.SHmax(:,j,k),'^','linewidth',1) %says that I cannot have more than two demensions for [S.SHmax] but seems to work well with Nmeasurements, 
      axis equal,xlim([-90,90]),ylim([-90,90]),grid
      xticks((-90:30:90)),yticks((-90:30:90))
      %legend('boreholes with breakouts','seismic stations with SWS fast directions',...
      %'location','northwest','FontSize',12,'FontWeight','bold')
      xlabel('SWS Fast Direction')
      ylabel('FMSHmax')
      %text(-88,-85,'d) ','FontSize',12,'FontWeight','bold')
      set(gcf,'renderer','Painters')
      title(['dmax = ',num2str(S.dmaxset(j)),'km zmax=',num2str(S.zmaxset(k)),'km'])
         end
      end
    
  %
  % Plot and save maps of Borehole SHmax vs FM SHmax estimates for various criteria
  %  - we did a version of this for SSA 2019, but haven't with the latest results
  %  - not expecting it to be particularly meaningful, but it seems to be helpful to
  %    have around, particularly for explaining the results in talks and such.
  %  - This takes ~25 seconds.  Once saved, no need to reproduce, just skip past.
  %
    sc=10; %how many km long should the SHmax stick be?

    figure(3),clf
    set(gcf,'Position',[100,400,665,440])
    %j=find(S.dmaxset==10); % this is the index for the 10km dmax solution, for testing
     for j=1:numel(S.dmaxset)

    [ha,pos]=tight_subplot(2,2,0.05,0.15,0.05); % this just helps it be a nicer looking journal figure
      axes(ha(1))
      % subplot(2,2,1),
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        % set(gca,'Color',[1 1 1]*0.75)
        hold on,
        plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',2)
        scatter(L.X,L.Y,50,L.FastDirection,'filled','markeredgecolor','k'),csym(90),colormap(cpolar),
        % colorbar('YTick',[-90:30:90])
        xticks([]),yticks([])
        % yticks(3730:20:3810)
        % title(['a) borehole SHmax'])
        text(227,3737,'a) SS FastDirection','FontSize',12,'FontWeight','bold')
        
      for k=1:numel(S.zmaxset) %using the 8 distances, it was set for zmaxset and still works for dmaxset
        %subplot(2,2,k+1)
        axes(ha(k+1)) % ran into a problem here, says exceeds the num of array elements(4)
        plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',L.X,L.Y,'ok'),axis equal,axis(P.R2)
        %set(gca,'Color',[1 1 1]*0.75)
        hold on,

        plot([L.X-sc*sind(S.SHmax(:,j,k)),L.X+sc*sind(S.SHmax(:,j,k))]',[L.Y-sc*cosd(S.SHmax(:,j,k)),L.Y+sc*cosd(S.SHmax(:,j,k))]','k','linewidth',2)%exceeds array bounds 
        scatter(L.X,L.Y,50,S.SHmax(:,j,k),'filled','markeredgecolor','k'),csym(90),colormap(cpolar),
        % colorbar('YTick',[-90:30:90])
        xticks([]),yticks([])
        % yticks(3730:20:3810)
        % title(['SHmax for all the [dmax = ',num2str(S.dmaxset(j)),'km] FM inversions: zmax=',num2str(S.zmaxset(k)),'km'])
        if k~=3
          text(227,3737,[char(97+k),') d_m_a_x = ',num2str(S.dmaxset(j)),' km, z_m_a_x = ',num2str(S.zmaxset(k)),' km'],'FontSize',12,'FontWeight','bold') % lowercase a = char(97)
        else
          text(227,3737,[char(97+k),') d_m_a_x = ',num2str(S.dmaxset(j)),' km, any depth'],'FontSize',12,'FontWeight','bold') % lowercase a = char(97)
        end
      end

      for k=1:4
        % subplot(2,2,k)
        axes(ha(k))
        h=image(RIutm(1:2),RIutm(3:4),flipud(I));
        uistack(h,'bottom')
      end
      h=colorbar('YTick',[-90:30:90],'Location','southoutside','FontSize',10);
      h.Position=[0.3 0.1 0.4 0.035];
      h.Label.String='SS FastDirection (deg E of N)';
      h.Label.FontSize=12;
      h.Label.FontWeight='bold';

       %figfilename=sprintf('figs/Mapfigures/map_FMSHmax_dmax%d',S.dmaxset(j));
       %saveas(gcf,figfilename,'epsc2')
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOW WELL DO BH FIT EACH OTHER IN THE GIVEN DISTANCE BINS?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Make a list of all possible pairs of BHs
  %  - calculate the distance and SHmax diffeence between boreholes
  %  - do this on a 2D grid (G), then reorder as as a list (L)
  %    Should be (57x57-57)/2 = 1596 unique pairs (full grid, minus diagonals, only one side)
  %
    %
    % Grid of locations of BH, grid of SHmax values, grid of SHmax uncertainties
    %
      [G.X1,G.X2]=meshgrid(B.X,B.X);
      [G.Y1,G.Y2]=meshgrid(B.Y,B.Y);
      [G.SHmax1,G.SHmax2]=meshgrid(B.SHmax,B.SHmax);
      [G.dSHmax1,G.dSHmax2]=meshgrid(B.dSHmax,B.dSHmax);
      [G.RegionNum1,G.RegionNum2]=meshgrid(L.RegionNum,L.RegionNum); %Regionnums were associated with SWS info
      
    %
    % Grid of index values, so we can identify the redundant ones, and convert
    % them to a list of pairs.  
    %
      [G.I,G.J]=meshgrid(1:numel(B.x),1:numel(B.x));
      i_redundant=find(G.I<=G.J);

      L.X1=G.X1(:);
      L.X2=G.X2(:);
      L.Y1=G.Y1(:);
      L.Y2=G.Y2(:);
      L.SHmax1=G.SHmax1(:);
      L.SHmax2=G.SHmax2(:);
      L.dSHmax1=G.dSHmax1(:);
      L.dSHmax2=G.dSHmax2(:);
      L.RegionNum1=G.RegionNum1(:);
      L.RegionNum2=G.RegionNum2(:);

      L.X1(i_redundant)=[];
      L.X2(i_redundant)=[];
      L.Y1(i_redundant)=[];
      L.Y2(i_redundant)=[];
      L.SHmax1(i_redundant)=[];
      L.SHmax2(i_redundant)=[];
      L.dSHmax1(i_redundant)=[];
      L.dSHmax2(i_redundant)=[];
      %L.RegionNum1(i_redundant)=[];
      %L.RegionNum2(i_redundant)=[];

    %
    % calculate the distance between the pairs, the difference in SHmax azimuth,
    % and the tolerance for their misfit (the sum of their two uncertainties)
    %
      L.dx=L.X1-L.X2;
      L.dy=L.Y1-L.Y2;
      L.d=sqrt(L.dx.^2+L.dy.^2); % distance in kilometers

      L.SHmaxDiff=daz(L.SHmax1,L.SHmax2);
    
    %
    % bin the BH v BH results, and calculate mean(d), rms, and N for each bin
    %  - use the same bins as the FM inversion dmax sets
    %  - actually, instead of doing bins, include everything up to that distance...
    %    this makes it more comparable to what the FM dmax means
    %
      for k=1:numel(S.dmaxset)
        i=find(L.d<=S.dmaxset(k));
        BvB.i{k,1}=i;
        BvB.N(k,1)=numel(i);
        BvB.d(k,1)=mean(L.d(i));
        BvB.rms(k,1)=rmsaz(L.SHmax1(i),L.SHmax2(i));
        BvB.macd(k,1)=mean(abs(L.SHmaxDiff(i))); % mean absolute circular difference
      end

    %
    % Can we get confidence intervals on the black line?  
    % Yes, by bootstrap resampling of BH SHmax values accounting for their uncertainty.
    % Overall this works, and gives results as expected, but when comparing to figure(9)a,
    % The PDFs are narrow enough that plotting them would confuse the journal figure 3a panel.
    % So, good that we did this, but really, we can leave it off.
    % 
    %
      Nbootstrap=1e4; % 1e3 is ok for binwidth 1, 1e5 for 0.2 binwidth for pub quality
      binwidth=.2;
      edges=binwidth/2:binwidth:90-binwidth/2;
      BvB.bins=binwidth:binwidth:90-binwidth;
      for k=3:numel(S.dmaxset)
        L.SHmax1draws=[rand(BvB.N(k),Nbootstrap).*L.dSHmax1(BvB.i{k})+L.SHmax1(BvB.i{k})];
        L.SHmax2draws=[rand(BvB.N(k),Nbootstrap).*L.dSHmax2(BvB.i{k})+L.SHmax2(BvB.i{k})];
        BvB.rmsPDF{k,1}=histcounts(rmsaz(L.SHmax1draws,L.SHmax2draws),edges,'normalization','pdf');
      end

      BvB.SCALE=2;
      figure(98),clf
      for k=3:numel(S.dmaxset)
        % plot(S.dmaxset(k)+BvB.rmsPDF{k},BvB.bins),hold on % this is just the lines, not the polygons
        polyy=[BvB.bins,BvB.bins(end),BvB.bins(1)];
        polyx=[S.dmaxset(k)+BvB.rmsPDF{k}*BvB.SCALE,S.dmaxset(k),S.dmaxset(k)];
        pgon=polyshape(polyx,polyy);
        plot(pgon,'edgecolor','k','facecolor',[1 1 1]*0.4)
        hold on
      end
      plot(S.dmaxset(3:end),BvB.rms(3:end),'k')
      plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
      ylim([0,90]),yticks(0:15:90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THE PDFs OF EACH BHvsFM COMPARISON, and some associated metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define binset and resolution on which to interpolate the FM and BH pdfs
  %
    binwidth=0.01; % this turns out to be the precision to which we can measure differences this way, so 2 sigfigs is good
    binwidth=0.1; % this one faster, with ok resolution.  Needs PDF SCALE adjust for plot
    binwidth=1; % this one is slightly faster, and bootstrap PDFs are a bit smoother
    interpbins=-90+binwidth:binwidth:90;

  %
  % Loop over each BH, each dmax, and each zmax
  % Calculate the difference PDFs:  This takes ~ 2 seconds, without bootstrapping
  %
    toc
    for i=1:numel(B.x)
    for j=1:numel(S.dmaxset)
    for k=1:numel(S.zmaxset)

  %
  % if FM solution doesn't exist, skip it and move on.  Actually, fill with NaNs (helps later)
  %
      if isnan(S.SHmax(i,j,k))
        BHvFMpdf{i,j,k}=interpbins+NaN;
        absBHvFMpdf{i,j,k}=interpbins+NaN;

        MU(i,j,k)=NaN;
        SIG(i,j,k)=NaN;
        P50(i,j,k)=NaN;
        P25(i,j,k)=NaN;
        P75(i,j,k)=NaN;
        P05(i,j,k)=NaN;
        P95(i,j,k)=NaN;
      else

  %
  % interpolate the two pdfs onto the same spacing, and re-Normalize them both
  %
        FMpdf=interp1(S.bins,S.SHmax_pdf{i,j,k},interpbins,'pchip','extrap');
        BHpdf=zeros(size(interpbins));
        BHpdf(find(interpbins>=B.SHmax(i)-B.dSHmax(i) & interpbins<=B.SHmax(i)+B.dSHmax(i)))=1;
        BHpdf(find(interpbins>=B.SHmax(i)-B.dSHmax(i)-180 & interpbins<=B.SHmax(i)+B.dSHmax(i)-180))=1;

        FMpdf=FMpdf/sum(FMpdf);
        BHpdf=BHpdf/sum(BHpdf);

  %
  % calculate the PDF of the difference using circular convolution
  % store the results, of both the wrapped difference and absolute difference
  %
        cDIFpdf=fftshift(cconv(BHpdf,fliplr(FMpdf),numel(interpbins)));
        cabsDIFpdf=cDIFpdf+fliplr(cDIFpdf);
        cabsDIFpdf(find(interpbins<=0))=0;

        BHvFMpdf{i,j,k}=cDIFpdf;
        absBHvFMpdf{i,j,k}=cabsDIFpdf;

  %
  % Make the metrics, while we're at it: mean, std, quartiles
  %
        MU(i,j,k)=sum(cabsDIFpdf.*interpbins);
        SIG(i,j,k)=sqrt(sum(abs(interpbins-sum(cabsDIFpdf.*interpbins)).^2.*cabsDIFpdf));

        cabsDIFcdf=cumsum(cabsDIFpdf); % cumulative distribution function
        [~,iquant]=min(abs(cabsDIFcdf-0.5));
        P50(i,j,k)=interpbins(iquant);

        [~,iquant]=min(abs(cabsDIFcdf-0.25));
        P25(i,j,k)=interpbins(iquant);
    
        [~,iquant]=min(abs(cabsDIFcdf-0.75));
        P75(i,j,k)=interpbins(iquant);

        [~,iquant]=min(abs(cabsDIFcdf-0.05));
        P05(i,j,k)=interpbins(iquant);
    
        [~,iquant]=min(abs(cabsDIFcdf-0.95));
        P95(i,j,k)=interpbins(iquant);

      end
    end
    end
    end
    toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE DIFFERENCE PDFs FOR INDIVIDUAL BOREHOLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % 6 panels total, 4 on LHS, 2 on RHS
  % 3 plots of seismogram-style wiggles of the individual PDFs (similar to previous Figure 2 or Figure 5)
  % 1 plot of FM depth distributions, along with CVM and breakout depths (New)
  % 1 plot of 5-95, 25-75, and 50 quantiles, as well as actual data (similar to previous Figure 6)
  % 1 Map (2D or 3D...) highlighting location of the individual borehole
  %
  % This takes 3.5 minutes for all of them, so only need to do it once to save the output
  %

  %
  % Plot the PDFs as seismo-wiggles: all 3 on the same subplot
  %
    % SCALE=3e3;
    SCALE=30/binwidth; % This one adapts for different bin widths
    warning('off','MATLAB:polyshape:repairedBySimplify') % keep polyshape from complaining every time....
    PCOLORS=[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250]*0.75; % these are the matlab-y ones
      
    i=randperm(57,1); % pick a random one to plot, for testing
    % for i=1:numel(B.x)
    figure(4),clf,
    set(gcf,'position',[840,140,550,820])

    subplot(312),
      for k=1:numel(S.zmaxset)
        for j=3:numel(S.dmaxset) % skip the "closest" ones
          %
          % skip it if there there's no FM solution
          %
          if ~isnan(S.SHmax(i,j,k))

            polyy=[interpbins(91:end),interpbins(end),interpbins(91)]; % just take the positive half...
            polyx=[S.dmaxset(j)+absBHvFMpdf{i,j,k}(91:end)*SCALE,S.dmaxset(j),S.dmaxset(j)];
            pgon=polyshape(polyx,polyy);
            plot(pgon,'edgecolor',PCOLORS(k,:),'linewidth',1,'facecolor',PCOLORS(k,:)*1.4)
            hold on
          end
        end
      end

      plot([0,40],[1 1]*90/sqrt(3),'k--','linewidth',1) % expected random rms
      ylim([0,90]),xlim([0,40])
      yticks(0:15:90)
      yticklabels({'0','','30','','60','','90'})
      xticks(0:5:35)
      set(gca,'FontSize',10)
      text(0.5,5,'b)','FontSize',12,'FontWeight','bold')
      xlabel('d_m_a_x (km)','FontSize',12,'FontWeight','bold'),ylabel('ACD (deg)','FontSize',12,'FontWeight','bold')

  %
  % depth distributions of FMs, as well as CVM depth and Breakout depths
  %  - oh, drat, don't have these as distributions, just the stats... oh well
  %
    subplot(313)
        polyx=[S.dmaxset(3:end)',fliplr(S.dmaxset(3:end)')];
        polyy=[S.z25th(i,3:end,k),fliplr(S.z75th(i,3:end,k))];
        
        polyx(find(isnan(polyy)))=[]; % if there are any NaNs, just skip over them
        polyy(find(isnan(polyy)))=[]; 

        pgon=polyshape(polyx,polyy);
        plot(pgon,'linestyle','none','facecolor',[1 1 1]*0.5,'HandleVisibility','off')
        hold on

        polyx=[S.dmaxset(3:end)',fliplr(S.dmaxset(3:end)')];
        polyy=[S.z05th(i,3:end,k),fliplr(S.z95th(i,3:end,k))];
        
        polyx(find(isnan(polyy)))=[]; % if there are any NaNs, just skip over them
        polyy(find(isnan(polyy)))=[]; 

        pgon=polyshape(polyx,polyy);
        plot(pgon,'linestyle','none','facecolor',[1 1 1]*0.75,'HandleVisibility','off')

        plot(S.dmaxset(3:end),S.z50th(i,3:end,3),'k','linewidth',1)

        axis ij

      plot([0,35],-[1 1]*B.z1el(i)/1e3,'r')
      plot([0,35],-[1 1]*B.z2el(i)/1e3,'r')
      pgon=polyshape([0,35,35,0],[-B.z1el(i),-B.z1el(i),-B.z2el(i),-B.z2el(i)]/1e3);
      plot(pgon,'linestyle','none','facecolor','r')

      %plot([0,35],-[1 1]*B.CVMz(i)/1e3,'k--','linewidth',1)

      xticks(0:5:35)
      set(gca,'FontSize',10)
      text(0.5,23,'c)','FontSize',12,'FontWeight','bold')


      ylim([-1,25]),xlim([0,40])
      xlabel('d_m_a_x (km)','FontSize',12,'FontWeight','bold'),
      ylabel('depth (km)','FontSize',12,'FontWeight','bold')

  %
  % Map showing at least where the borehole is
  %
    subplot(311)
      plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w',B.X,B.Y,'ok'),axis equal,axis(P.R2)
      hold on,
      plot([B.X(i)-sc*sind(B.SHmax(i)),B.X(i)+sc*sind(B.SHmax(i))]',[B.Y(i)-sc*cosd(B.SHmax(i)),B.Y(i)+sc*cosd(B.SHmax(i))]','k','linewidth',2)
      scatter(B.X(i),B.Y(i),50,B.SHmax(i),'filled','markeredgecolor','k'),colorbar,csym(90),colormap(cpolar),
      h=colorbar('YTick',[-90:30:90],'FontSize',10);
      h.Label.String='SHmax (deg E of N)';
      h.Label.FontSize=12;
      h.Label.FontWeight='bold';

      yticks([]),xticks([])
  
      h=image(RIutm(1:2),RIutm(3:4),flipud(I));
      uistack(h,'bottom')
      text(227,3737,['a) borehole #',num2str(i)],'FontSize',12,'FontWeight','bold')

  %
  % Save the figure, per borehole
  %
    % figfilename=sprintf('figs/supplementalfigs/BH%02dsummary',i);
    % saveas(gcf,figfilename,'epsc2')
    % end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STACK EACH DMAX and ZMAX OVER THE ENTIRE REGION (RMS PLOT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Calculate the expected PDF of the sum of BHs (so, when normalized, the mean), using the 
  %   convolution theorem in the Fourier domain. Loop over all the zmax and dmax criteria
  %
    bins=interpbins(find(interpbins>0)); % these are the right-half bins, for absolute difference [0,90]

  %
  % Draw Lots of samples from each difference PDF to be stacked, using ranany
  %  Then square, sum, divide by relevant N, and sqrt to get RMS
  % Can this be fast?
  %
      Nbootstrap=1e4; % how many samples to draw
      edges=[bins-binwidth/2,bins(end)+binwidth/2]; % for use in estimating PDF

      toc
      for k=1:numel(S.zmaxset)
      for j=1:numel(S.dmaxset)

      %
      % identify the non-nan solutions for this set of dmax and zmax
      % Populate an array NBHxNbootstrap from which to calculate the RMS
      % Calculate NBH realizations of RMS
      %
          NBH=sum(~isnan(S.SHmax(:,j,k))); % This is the number of BH for which FM solutions exist, [0,57]
          iBH=find(~isnan(S.SHmax(:,j,k))); % This is the index of BH for which FM solutions exist

        %
        % if FM solution doesn't exist, skip it and move on.  Actually, fill with NaNs (helps later)
        %
          if NBH==0
            bootRMSdraws{j,k}=zeros(1,Nbootstrap)*NaN;
            bootRMSpdf{j,k}=zeros(1,numel(bins))*NaN;
          else

            clear bootdraws
            for i=1:NBH
              bootdraws(i,:)=ranany(interpbins,absBHvFMpdf{iBH(i),j,k},Nbootstrap)';
            end
        
            bootRMSdraws{j,k}=sqrt(sum(bootdraws.^2,1)/NBH);
            bootRMSpdf{j,k}=histcounts(bootRMSdraws{j,k},edges,'normalization','probability');
          end
        end
        toc
        end

  %
  % Calculate the meaningful metrics to be plotted
  %
    bMU=cell2mat(cellfun(@mean,bootRMSdraws,'uniformoutput',false));
    bSIG=cell2mat(cellfun(@std,bootRMSdraws,'uniformoutput',false));

    bootRMScdf=cellfun(@cumsum,bootRMSpdf,'uniformoutput',false); % cumulative distribution function

    [~,iquant]=cellfun(@(x) min(abs(x-0.05)),bootRMScdf);
    bP05=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.25)),bootRMScdf);
    bP25=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.50)),bootRMScdf);
    bP50=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.75)),bootRMScdf);
    bP75=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.95)),bootRMScdf);
    bP95=bins(iquant);

    %Re-NaN the points for solutions that don't exist (otherwise, they're zero)
    bP05(find(isnan(bMU)))=NaN; 
    bP25(find(isnan(bMU)))=NaN; 
    bP50(find(isnan(bMU)))=NaN; 
    bP75(find(isnan(bMU)))=NaN; 
    bP95(find(isnan(bMU)))=NaN; 

    % RMS over all the BHvFM differences
    BHvFMrms=sqrt(reshape(sum(BHvFM.^2,'omitnan'),16,3)'./reshape(sum(~isnan(BHvFM)),16,3)')';

  %
  % Plot the results with reconstituted RMS
  %
    %
    % Plot the PDFs as seismo-wiggles individually, then quantile swaths below
    %
      SCALE=10/binwidth;
      figure(5),clf,
        for k=1:numel(S.zmaxset)
          subplot(4,1,k),
          for j=3:numel(S.dmaxset) % skip the "closest" ones
            %
            % skip it if there there's no FM solutions for this criteria
            %
            if sum(~isnan(S.SHmax(:,j,k)))>0

              polyy=[bins,bins(end),bins(1)];
              polyx=[S.dmaxset(j)+bootRMSpdf{j,k}*SCALE,S.dmaxset(j),S.dmaxset(j)];
              pgon=polyshape(polyx,polyy);
              plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))
              hold on
            end
          end

          plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
          plot(S.dmaxset(3:end),BHvFMrms(3:end,k),'k.-','linewidth',1) % the RMS instead of the mean

          ylim([25,55]),%yticks(0:30:90)
          xlim([0,40]),grid,
          title(sprintf('Stacked over entire region, zmax = %d km',S.zmaxset(k)))
          xlabel('dmax (km)'),ylabel('RMS (deg)')
        end

        subplot(414)
        for k=1:numel(S.zmaxset)
          polyx=[S.dmaxset(3:end)',fliplr(S.dmaxset(3:end)')];
          polyy=[bP25(3:end,k);flipud(bP75(3:end,k))]';
          
          polyx(find(isnan(polyy)))=[]; % if there are any NaNs, just skip over them
          polyy(find(isnan(polyy)))=[]; 

          pgon=polyshape(polyx,polyy);
          plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))
          hold on

          polyx=[S.dmaxset(3:end)',fliplr(S.dmaxset(3:end)')];
          polyy=[bP05(3:end,k);flipud(bP95(3:end,k))]';
          
          polyx(find(isnan(polyy)))=[]; % if there are any NaNs, just skip over them
          polyy(find(isnan(polyy)))=[]; 

          pgon=polyshape(polyx,polyy);
          plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))

          % plot(S.dmaxset(3:end),bP50(3:end,k),'color',PCOLORS(k,:)/0.75)
        end

        set(gca,'colororderindex',1)
        plot(S.dmaxset(3:end),bP50(3:end,:))
        plot(S.dmaxset(3:end),BvB.rms(3:end),'k:','linewidth',1) % How well do BH fit one another across entire region?

        plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
        plot(xlim,[1 1]*rmsaz(B.SHmax,B.FMSHmax),'k:') % the YH13 v BH rms

        ylim([25,55]),%yticks(0:15:90)
        xlim([0,35]),grid
        title('PDF of RMS over entire region')
        xlabel('dmax (km)'),ylabel('RMS (deg)')

    %
    % Plot the PDFs as seismo-wiggles stacked on the same plot
    %
      figure(6),clf,
        SCALE=10/binwidth*0.5;
        for k=1:numel(S.zmaxset)
          for j=3:numel(S.dmaxset) % skip the "closest" ones
            %
            % skip it if there there's no FM solutions for this criteria
            %
            if sum(~isnan(S.SHmax(:,j,k)))>0

              % SCALE=2/max(bootRMSpdf{j,k}); % whatever it takes to get the pdf to be 2km high...
              polyy=[bins,bins(end),bins(1)];
              polyx=[S.dmaxset(j)+bootRMSpdf{j,k}*SCALE,S.dmaxset(j),S.dmaxset(j)];
              pgon=polyshape(polyx,polyy);
              % plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))
              plot(pgon,'edgecolor',PCOLORS(k,:),'facecolor',PCOLORS(k,:))
              % plot(pgon,'edgecolor',PCOLORS(k,:),'facecolor','none')
              hold on
            end
          end

          plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
          % plot(S.dmaxset(3:end),BHvFMrms(3:end,k),'-','color',PCOLORS(k,:)) % the RMS instead of the mean
          % plot(S.dmaxset(3:end),bP50(3:end,k),'-','color',PCOLORS(k,:)) % the median values

          ylim([15,75]),yticks(0:15:90)
          xlim([0,40]),grid,
          % title(sprintf('Stacked over entire region, zmax = %d km',S.zmaxset(k)))
          xlabel('dmax (km)'),ylabel('RMS (deg)')
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SUBREGIONS BY GEOLOGIC CONTEXT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Previously, this was calculated here.  Now, it's just being loaded from TableS1 as part of
  % load_and_parse_BHFM.m subroutine.
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
%       this info excludes the area that isnt the LAbasin   
  
  %
  % plot to check and make sense out of all of this
  %
    LINES=lines(5);

    figure(7),clf, %editing to add SS in 2D
    scatter(F.X,F.Y,5,-F.z,'filled'),
    hold on
    plot(P.cutm(:,1),P.cutm(:,2),'k'),
    plot(P.futm(:,1),P.futm(:,2),'color',[1 1 1]*0.5),axis equal,axis(P.R2),
    contour(P.CVM.X/1e3,P.CVM.Y/1e3,P.CVM.Z/1e3,[-20:2:1],'linewidth',1),colorbar,caxis([-20,1])
    plot(F.X(iNorthRidgeEQ),F.Y(iNorthRidgeEQ),'k*','markersize',20,'linewidth',2) %FM
    plot([B.X-sc*sind(B.SHmax),B.X+sc*sind(B.SHmax)]',[B.Y-sc*cosd(B.SHmax),B.Y+sc*cosd(B.SHmax)]','k','linewidth',2) %BH
    %plot([L.X(iLA)-sc*sind(L.FastDirection(iLA)),L.X(iLA)+sc*sind(L.FastDirection(iLA))]',[L.Y(iLA)-sc*cosd(L.FastDirection(iLA)),L.Y(iLA)+sc*cosd(L.FastDirection(iLA))]','k','linewidth',0.05) %SS
    %scatter(L.X(iLA),L.Y(iLA),15,L.ResultantLength(iLA),'filled'),colorbar,grid %csym(90)
    plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',0.05) %SS
    scatter(L.X,L.Y,10,L.ResultantLength,'filled'),colorbar,grid %csym(90)
    for i=1:max(B.GeolNum)
      plot(B.X(find(B.GeolNum==i)),B.Y(find(B.GeolNum==i)),'o','markeredgecolor','k','markerfacecolor',LINES(i,:))
    end
    colormap(jet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMS PLOT FOR GEOLOGIC CONTEXT SUBREGIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % First plot the map with BH colored by their geologic context.
  % Actually, for figure 9, try out having the pdfs of the entire region in one figure.  
  %
    figure(8),clf
    subplot(3,2,1)
      plot(P.cutm(:,1),P.cutm(:,2),'k',P.futm(:,1),P.futm(:,2),'w'),axis equal,axis(P.R2)
      set(gca,'Color',[1 1 1]*0.75)
      hold on,
      % plot([B.X(iB)-sc*sind(B.SHmax(iB)),B.X(iB)+sc*sind(B.SHmax(iB))]',[B.Y(iB)-sc*cosd(B.SHmax(iB)),B.Y(iB)+sc*cosd(B.SHmax(iB))]','k','linewidth',2)
      scatter(B.X,B.Y,50,B.GeolNum,'filled'),colorbar,colormap(lines(5)),caxis([0.5,5.5])
      plot(B.X,B.Y,'ok')
      title(['Geologic Context'])

    figure(9),clf
    subplot(3,2,1)
        SCALE=10/binwidth;
        for k=1:numel(S.zmaxset)
          for j=3:numel(S.dmaxset) % skip the "closest" ones
            %
            % skip it if there there's no FM solutions for this criteria
            %
            if sum(~isnan(S.SHmax(:,j,k)))>0

              % SCALE=2/max(bootRMSpdf{j,k}); % whatever it takes to get the pdf to be 2km high...
              polyy=[bins,bins(end),bins(1)];
              polyx=[S.dmaxset(j)+bootRMSpdf{j,k}*SCALE,S.dmaxset(j),S.dmaxset(j)];
              pgon=polyshape(polyx,polyy);
              % plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))
              plot(pgon,'edgecolor',PCOLORS(k,:),'facecolor',PCOLORS(k,:)*1.4)
              % plot(pgon,'edgecolor',PCOLORS(k,:),'facecolor','none')
              hold on
            end
          end

          plot(S.dmaxset(3:end),BHvFMrms(3:end,k),'-','color',PCOLORS(k,:)) % the rms of the best fitting values
          % plot(S.dmaxset(3:end),bP50(3:end,k),'-','color',PCOLORS(k,:)) % the median values
        end

        plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
        plot(S.dmaxset(3:end),BvB.rms(3:end),'k') % How well do BH fit one another across entire region?

        ylim([0,90]),yticks(0:15:90)
        xlim([0,40]),%grid,
        title(sprintf('entire region, N=57'))
        xlabel('dmax (km)'),ylabel('RMS (deg)')

      %
      % Optionally Add the PDFs of the BHvBH line.  These are small enough to be not that helpful in this plot,
      % so leave them off to avoid confusing the issue.
      %
        for j=3:numel(S.dmaxset)
          polyy=[BvB.bins,BvB.bins(end),BvB.bins(1)];
          polyx=[S.dmaxset(j)+BvB.rmsPDF{j}*BvB.SCALE,S.dmaxset(j),S.dmaxset(j)];
          pgon=polyshape(polyx,polyy);
          h=plot(pgon,'edgecolor','k','facecolor',[1 1 1]*0.4);
          uistack(h,'bottom')
        end

  
  %
  % Loop through all subregions, plot their results on a single figure
  %
    % iR=find(strcmp(Geol.name,'Basin Edge'));
    for iR=1:max(B.GeolNum)
    iB=find(B.GeolNum==iR);

  %
  % Draw Lots of samples from each difference PDF to be stacked, using ranany
  %  Then square, sum, divide by relevant N, and sqrt to get RMS
  %
    Nbootstrap=1e4; % how many samples to draw: 1e4 is ok, 1e5 smooths out the lumps, 1e6 moreso but takes 2 minutes
    edges=[bins-binwidth/2,bins(end)+binwidth/2]; % for use in estimating PDF

    toc
    for k=1:numel(S.zmaxset)
    for j=1:numel(S.dmaxset)

    %
    % identify the non-nan solutions for this set of dmax and zmax
    % Populate an array NBHxNbootstrap from which to calculate the RMS
    % Calculate NBH realizations of RMS
    %
      NBH=sum(~isnan(S.SHmax(iB,j,k))); % This is the number of BH for which FM solutions exist, [0,57]
      iBH=find(~isnan(S.SHmax(iB,j,k))); % This is the index of BH for which FM solutions exist

    %
    % if FM solution doesn't exist, skip it and move on.  Actually, fill with NaNs (helps later)
    %
      if NBH==0
        sbootRMSdraws{j,k}=zeros(1,Nbootstrap)*NaN;
        sbootRMSpdf{j,k}=zeros(1,numel(bins))*NaN;
      else

        clear bootdraws
        for i=1:NBH
          bootdraws(i,:)=ranany(interpbins,absBHvFMpdf{iB(iBH(i)),j,k},Nbootstrap)';
        end
    
        sbootRMSdraws{j,k}=sqrt(sum(bootdraws.^2,1)/NBH);
        sbootRMSpdf{j,k}=histcounts(sbootRMSdraws{j,k},edges,'normalization','probability');
      end

    end
    toc
    end

  %
  % Calculate the meaningful metrics to be plotted
  %
    sbMU=cell2mat(cellfun(@mean,sbootRMSdraws,'uniformoutput',false));
    sbSIG=cell2mat(cellfun(@std,sbootRMSdraws,'uniformoutput',false));

    sbootRMScdf=cellfun(@cumsum,sbootRMSpdf,'uniformoutput',false); % cumulative distribution function

    [~,iquant]=cellfun(@(x) min(abs(x-0.05)),sbootRMScdf);
    sbP05=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.25)),sbootRMScdf);
    sbP25=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.50)),sbootRMScdf);
    sbP50=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.75)),sbootRMScdf);
    sbP75=bins(iquant);

    [~,iquant]=cellfun(@(x) min(abs(x-0.95)),sbootRMScdf);
    sbP95=bins(iquant);

    %Re-NaN the points for solutions that don't exist (otherwise, they're zero)
    sbP05(find(isnan(sbMU)))=NaN; 
    sbP25(find(isnan(sbMU)))=NaN; 
    sbP50(find(isnan(sbMU)))=NaN; 
    sbP75(find(isnan(sbMU)))=NaN; 
    sbP95(find(isnan(sbMU)))=NaN; 

    % RMS over all the BHvFM differences
    sBHvFMrms=sqrt(reshape(sum(BHvFM(iB,:,:).^2,1,'omitnan'),16,3)'./reshape(sum(~isnan(BHvFM(iB,:,:)),1),16,3)')';


  %
  % Plot the results with reconstituted RMS, parallel to Figure 5 above
  %
    %
    % Plot the RMS quantile vs dmax shapes, along with "best" data points
    %
      figure(8)
      subplot(3,2,iR+1)
        for k=1:numel(S.zmaxset)
          polyx=[S.dmaxset(3:end)',fliplr(S.dmaxset(3:end)')];
          polyy=[sbP25(3:end,k);flipud(sbP75(3:end,k))]';
          
          polyx(find(isnan(polyy)))=[]; % if there are any NaNs, just skip over them
          polyy(find(isnan(polyy)))=[]; 

          pgon=polyshape(polyx,polyy);
          plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))
          hold on

          polyx=[S.dmaxset(3:end)',fliplr(S.dmaxset(3:end)')];
          polyy=[sbP05(3:end,k);flipud(sbP95(3:end,k))]';
          
          polyx(find(isnan(polyy)))=[]; % if there are any NaNs, just skip over them
          polyy(find(isnan(polyy)))=[]; 

          pgon=polyshape(polyx,polyy);
          plot(pgon,'linestyle','none','facecolor',PCOLORS(k,:))

          % plot(S.dmaxset(3:end),sbP50(3:end,k),'color',PCOLORS(k,:)/0.75)
        end

        set(gca,'colororderindex',1)
        plot(S.dmaxset(3:end),sbP50(3:end,:))
        plot(S.dmaxset(3:end),BvB.rms(3:end),'k:') % How well do BH fit one another across entire region?

        plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
        plot(xlim,[1 1]*rmsaz(B.SHmax,B.FMSHmax),'k:') % the YH13 v BH rms

        ylim([0,90]),yticks(0:15:90)
        xlim([0,35]),grid
        title(sprintf('%s',B.Geol{iB(1)}))
        xlabel('dmax (km)'),ylabel('RMS (deg)')

    %
    % Plot the PDFs as seismo-wiggles, but instead of a map in the first place, include the entire region results.
    %
      figure(9)
      subplot(3,2,iR+1),
      SCALE=10/binwidth;
        for k=1:numel(S.zmaxset)
          for j=3:numel(S.dmaxset) % skip the "closest" ones
            %
            % skip it if there there's no FM solutions for this criteria
            %
            if sum(~isnan(S.SHmax(iB,j,k)))>0

              polyy=[bins,bins(end),bins(1)];
              polyx=[S.dmaxset(j)+sbootRMSpdf{j,k}*SCALE,S.dmaxset(j),S.dmaxset(j)];
              pgon=polyshape(polyx,polyy);
              plot(pgon,'edgecolor',PCOLORS(k,:),'facecolor',PCOLORS(k,:)*1.4)
              hold on
            end
          end

          plot(S.dmaxset(3:end),sBHvFMrms(3:end,k),'color',PCOLORS(k,:)) % RMS of best estimates for BH in this Geol Context
        end
        plot(xlim,[1 1]*90/sqrt(3),'k--') % expected random rms
        ylim([0,90]),yticks(0:15:90)
        xlim([0,40]),%grid,
        title(sprintf('%s, N=%d',B.Geol{iB(1)},numel(iB)))
        xlabel('dmax (km)'),ylabel('RMS (deg)')
      set(gcf,'renderer','Painters')


  %
  % Save the figure, per subregions
  %
    % figfilename=sprintf('figs/map_FMSHmax_values/SubRegionGeol%sSummaryRMS',Geol.shortname{iR});
    % set(gcf,'InvertHardcopy','off')
    % saveas(gcf,figfilename,'epsc2')
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE ILLUSTRATING DEPTH DEPENDENCE WITH BH # 52-55 (now 21,22,52,23) AT NEWPORT-INGLEWOOD FAULT WLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Load/create extra stuff we'll need
  %  - CFM fault traces (different from USGS historic ruptures)
  %  - CFM 3D surface for Newport-Inglewood fault
  %  - which BH are we focused on, and what extent (+/- 5 km from the BH)
  %  - Define a profile across and sample the CVM and topo at those locations for plotable profiles
  %  - borehole breakout midpoint, for picking a color...
  %  - figure out the 5/95 and 25/75 pdf range for the FMs.  
  %    Use zmax = any:  there aren't any shallow earthquakes within 5+ km
  %    Use dmax = 4 km: the one with solutions existing for all 4 boreholes for any depth
  %    52-54 are really really the same, 55 is very close... So, just use the 52 one
  %
    %
    % CFM stuff for parts 1 and 2
    %
      % C=LoadCFMxyznisidi_v5_2('*NIRC-LABS*'); % Load newport-inglewood fault in 3D
      % P.NIF.x=C{2}(:,1);
      % P.NIF.y=C{2}(:,2);
      % P.NIF.z=C{2}(:,3);
      % [P.NIF.X,P.NIF.Y]=ll2utm(P.NIF.y,P.NIF.x); 

      % filename='/Users/kluttrell/CaliforniaStress/CommunityModels/CFM5_release_2017/obj/CFM5_preferred_traces/CFM5_preferred_traces.lonLat';
      % [~,r]=system(['cat ',filename,' | grep -v \# | sed ''s/^>.*/NaN NaN NaN/''']);
      % C=textscan(r,'%f %f %f');
      % P.CFM.x=C{1};
      % P.CFM.y=C{2};
      % P.CFM.z=C{3};
      % [P.CFM.X,P.CFM.Y]=ll2utm(P.CFM.y,P.CFM.x); 
    
    %
    % Define the boreholes and region of interest
    % 
      i=[21,22,52,23];
      P.R4=[372,382.1,3754,3763.5];
      sqrt(sum([B.X(i)-B.X(i(1)),B.Y(i)-B.Y(i(1))].^2,2)) % this is how far apart they are

    %
    % Crosssection Profile stuff for part 2
    %
      profileA.ends=[372.1,3755.5;382.1,3761.2];
      profileA.X=[linspace(profileA.ends(1,1),profileA.ends(2,1),1e2)'];
      profileA.Y=[linspace(profileA.ends(1,2),profileA.ends(2,2),1e2)'];
      [profileA.y,profileA.x]=utm2ll(profileA.X*1e3,profileA.Y*1e3,11);
      profileA.d=[0;cumsum(sqrt(sum(diff([profileA.X,profileA.Y]).^2,2)))];
      
      fid=fopen('temp','w');
      fprintf(fid,'%f %f\n',[profileA.X,profileA.Y]'*1e3);
      fclose(fid);
      [~,r]=system('cat temp | gmt grdtrack -Ggrds/CVM15_250m_UTM11.grd -Z');
      profileA.CVMz=cell2mat(textscan(r,'%f'));

      fid=fopen('temp','w');
      fprintf(fid,'%f %f\n',[profileA.x,profileA.y]');
      fclose(fid);
      [~,r]=system('cat temp | gmt grdtrack -Ggrds/topo.grd -Z');
      profileA.topoz=cell2mat(textscan(r,'%f'));
      system('rm temp');

      profileA.azimuth=atand((profileA.ends(2,1)-profileA.ends(1,1))/(profileA.ends(2,2)-profileA.ends(1,2)));

      Bdi=sqrt(sum([B.X(i)-profileA.ends(1,1),B.Y(i)-profileA.ends(1,2)].^2,2));

      B.zmidel=mean([B.z1el,B.z2el],2);

    %
    % Local FMs rotated for Crosssection View plot in part 2
    %  - find the subset of FMs, so we don't have to rotate them all
    %  - define the rotation matrix based on the azimuth calculated above
    %
      iFM=find(F.X>P.R4(1) & F.X<P.R4(2) & F.Y>P.R4(3) & F.Y<P.R4(4)); % the ones in range

      theta=90-profileA.azimuth; % rotation amount
      Rot=[cosd(theta),sind(theta);-sind(theta),cosd(theta)];

      FXrot=(F.X(iFM)-profileA.ends(1,1))*Rot(1,1)+(F.Y(iFM)-profileA.ends(1,2))*Rot(1,2);
      FYrot=(F.X(iFM)-profileA.ends(1,1))*Rot(2,1)+(F.Y(iFM)-profileA.ends(1,2))*Rot(2,2);
      FZrot=-F.z(iFM);

    %
    % FM quantiles for part 3
    %
      FMpdf=interp1(S.bins,S.SHmax_pdf{i(1),find(S.dmaxset==6),find(S.zmaxset==999)},interpbins,'pchip','extrap');
      FMcdf=cumsum(FMpdf);
      FMcdf=FMcdf/max(FMcdf);

      [~,iquant]=min(abs(FMcdf-0.05));
      FMp05=interpbins(iquant);

      [~,iquant]=min(abs(FMcdf-0.25));
      FMp25=interpbins(iquant);

      [~,iquant]=min(abs(FMcdf-0.50));
      FMp50=interpbins(iquant);

      [~,iquant]=min(abs(FMcdf-0.75));
      FMp75=interpbins(iquant);

      [~,iquant]=min(abs(FMcdf-0.95));
      FMp95=interpbins(iquant);


  %
  % Make a 3D one, just so we can make sure we understand what's going on...
  %
    figure(10),clf,
    scatter3(F.X,F.Y,-F.z,10,-F.z,'filled'),colorbar,caxis([-20,0])
    axis equal,axis(P.R4),zlim([-20,1])
    hold on
    plot(B.X(i),B.Y(i),'ko','linewidth',1)
    plot3((B.X(i)*[1 1])',(B.Y(i)*[1 1])',[B.z1el(i) B.z2el(i)]'/1e3,'-r','linewidth',3)
    surf(P.CVM.X(1:2:end)/1e3,P.CVM.Y(1:2:end)/1e3,P.CVM.Z(1:2:end,1:2:end)/1e3,'facecolor','none') 
    % plot3(P.CFM.X/1e3,P.CFM.Y/1e3,P.CFM.z/1e3,'k')
    % plot3(P.NIF.X/1e3,P.NIF.Y/1e3,P.NIF.z/1e3,'k.')
    % plot([L.X-sc*sind(L.FastDirection),L.X+sc*sind(L.FastDirection)]',[L.Y-sc*cosd(L.FastDirection),L.Y+sc*cosd(L.FastDirection)]','k','linewidth',0.05) %SS
    % scatter(L.X,L.Y,15,L.ResultantLength,'filled'),colorbar,grid %csym(90)
    contour(P.CVM.X/1e3,P.CVM.Y/1e3,P.CVM.Z/1e3,1e2,'linewidth',0.5),
    plot(profileA.ends(:,1),profileA.ends(:,2),'r')
    plot3(profileA.X,profileA.Y,profileA.topoz/1e3,'r')
    plot3(profileA.X,profileA.Y,profileA.CVMz/1e3,'r')

    % view(90-profileA.azimuth,0) % cross-section view

  %
  % Figure Part 1: mapview and context
  %
    COLORS=lines(4);
    
    figure(11),clf
    subplot(131),

    scatter(F.X,F.Y,5,-F.z,'filled'),colorbar,caxis([-20,0])
    axis equal,axis(P.R4)
    hold on,
    contour(P.CVM.X/1e3,P.CVM.Y/1e3,P.CVM.Z/1e3,[-10:0.5:0],'linewidth',0.5),
    
    % plot(B.X(i),B.Y(i),'ko','linewidth',1)
    for k=1:numel(i)
      scatter(B.X(i(k)),B.Y(i(k)),50,COLORS(k,:),'filled')
    end

    % plot(P.CFM.X/1e3,P.CFM.Y/1e3,'k')
    % plot(P.NIF.X/1e3,P.NIF.Y/1e3,'k.')
    plot(profileA.ends(:,1),profileA.ends(:,2),'r')

  %
  % Figure Part 2: cross section view and context
  %

    subplot(132)
    plot(profileA.d,profileA.topoz/1e3,'k',profileA.d,profileA.CVMz/1e3,'k')
    axis equal,axis tight
    xlabel('distance along profile (km)'),ylabel('depth')
    hold on,
      % plot(Bdi,B.topoz(i)/1e3,'ko','linewidth',1)
      % plot((Bdi*[1 1])',[B.z1el(i) B.z2el(i)]'/1e3,'-r','linewidth',3)
    for k=1:numel(i)
      plot((Bdi(k)*[1 1])',[B.z1el(i(k)) B.z2el(i(k))]'/1e3,'color',COLORS(k,:),'linewidth',3)
    end

    halfwidth=0.5; % how wide do you want the fault swath to be, relative to the boreholes?
    polyx=mean(Bdi(1:3))+[-1 1 1 -1]*halfwidth;
    [~,iNIFd1]=min(abs(profileA.d-(mean(Bdi(1:3))-halfwidth)));
    [~,iNIFd2]=min(abs(profileA.d-(mean(Bdi(1:3))+halfwidth)));
    polyy=[profileA.topoz([iNIFd1,iNIFd2]);profileA.CVMz([iNIFd2,iNIFd1])]'/1e3;
    pgon=polyshape(polyx,polyy);
    plot(pgon,'linestyle','none','facecolor',[1 1 1]*0.75)

    scatter(FXrot,FZrot,10,FZrot,'filled')
    colorbar,caxis([-20,0])

  %
  % Figure Part 3: azimuth illustration
  %

    subplot(133),
    for k=1:numel(i)
      polarplot(d2r(B.SHmax(i(k)))*[1 1],[-1,1],'color',COLORS(k,:),'linewidth',3)
      polarplot(d2r(B.SHmax(i(k))+[-1 1 1 -1 -1]*B.dSHmax(i(k))),[1 1 -1 -1 1],'color',COLORS(k,:)*0.75)
      hold on
    end
    set(gca,'ThetaDir','clockwise','ThetaZeroLocation','top','Rtick',0)

    polarplot(d2r(S.SHmax(i(1),find(S.dmaxset==6),find(S.zmaxset==999)))*[1 1],[-1,1],'k','linewidth',1)
    polarplot(d2r([FMp05 FMp95 FMp95 FMp05 FMp05]),[1 1 -1 -1 1],'k')
    % polarplot(d2r([FMp25 FMp75 FMp75 FMp75 FMp25]),[1 1 -1 -1 1],'k')
    % title('SHmax (degEofN)')

    set(gcf,'renderer','Painters')


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT TO ANSWER THE NORTHRIDGE ABUNDACE QUESTION, because I keep making it anyways...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % which quakes are within the region? (well, we need which ones are within 35 km, but oh well)
  %
    i=find(F.x>P.R(1)&F.x<P.R(2)&F.y>P.R(3)&F.y<P.R(4));
    figure(99),clf
    subplot(211),plot(F.t,F.M,'ok',F.t(i),F.M(i),'.r'),datetick
    legend('all quakes','within region','location','northwest'),ylabel('M')
    subplot(212),histogram(F.t),hold on,histogram(F.t(i)),datetick
    legend('all quakes','within region','location','northwest'),ylabel('N')



toc
