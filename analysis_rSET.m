clear all;
%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load in Louisiana coast map data
% load the CRMS data quickly from the mat file
% (see CRMS_Surface_Elevation_Internship.m for *.csv to *.mat)
%
	
	filename=['Script&Data/LAparishes.xy'];
	fid3=fopen(filename);
	C=textscan(fid3,'%f %f');
	fclose(fid3);

		parishlon=C{1};
		parishlat=C{2};

	filename=['Script&Data/coast.xy'];
	fid4=fopen(filename);
	C=textscan(fid4,'%f %f');
	fclose(fid4);

		coastlon=C{1};
		coastlat=C{2};

	load CRMS_rSET_heights.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	
% plotting station lat/longs to visualize locations spatially 
%	

	figure(1),
	hold on;
	plot(coastlon,coastlat,'k',parishlon,parishlat,'k');
	s1=scatter(stationlon,stationlat,'.b');
	set(gca,'DataAspectRatio',[1 1 1/cosd(30)])
	xlabel('Longitude')
	ylabel('Latitude')
	title('Station Map')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plotting rSET mean pin height for all CRMS stations
%

	Hrel=cellfun(@(x) x-x(1),H,'uniformoutput',false);

	%plotting absolute mean pin heights
%{
	figure(2),clf
	plot(cell2mat(T),cell2mat(H),'*'); datetick;
	ylabel('Absolute rSET Mean Height (mm)')
	xlabel('Time')
	title('Absolute Mean rSET Height for all CRMS Stations (Blue) vs. South Stations (Orange)')

%}

	
	%what about just the south ones?
	iSouth=find(stationlat<29.5);
	hold on;
	plot(cell2mat(T(iSouth)),cell2mat(H(iSouth)),'.')
	hold off;

	%{
%plotting relative mean pin heights
	figure(3),clf;
	plot(cell2mat(T),cell2mat(Hrel),'*'); datetick;
	ylabel('Mean rSET Height (mm), Relative to First Measurement')
	xlabel('Time')
	title('Relative Mean rSET Height for all CRMS Stations (Blue) vs. South Stations (Orange)')
	hold on; 
	plot(cell2mat(T(iSouth)),cell2mat(Hrel(iSouth)),'.')
	hold off;
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plotting individual CRMS station data
%

%{
	for i=1:numel(stationset);
		
		figure(4), clf;
		plot(T{i},Hrel{i},'*'); datetick;
		xlabel('Time')
		ylabel('Relative Mean Pin Height (mm)')
		title(stationset{i});
		grid on;

		pause

	end 
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Loading in GPS Station Data and plotting on map
%

		%loading in GPS lat/longs	
		filename=['GPS_All/gpslocations.txt'];
		fid2=fopen(filename);
		C=textscan(fid2,'%s %f %f');
		fclose(fid2);
		
			GPSnames=C{1};		
			GPSlat=C{2};
			GPSlong=C{3};

			nGPS=length(GPSlat);

	for i=1:nGPS;

		filename=['GPS_All/',GPSnames{i},'.txt'];
		fid=fopen(filename);
		C=textscan(fid,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',1);
		fclose(fid);

		%
		%assigning variables to the columns of data
		%

		rawGPS.stationid{i,1}=C{1}(1);
		rawGPS.t{i,1}=datenum(C{2},'yymmmdd');
		rawGPS.tyear{i,1}=C{3};
		%rawGPS.x{i,1}=(C{8}+C{9});
		%rawGPS.y{i,1}=(C{10}+C{11});
		rawGPS.z{i,1}=(C{12}+C{13});
		%rawGPS.sigx{i,1}=C{15};
		%rawGPS.sigy{i,1}=C{16};
		rawGPS.sigz{i,1}=C{17};

		%
		%zero data to first of time series(changing units to mm)
		%

		%rawGPS.x{i,1}=(rawGPS.x{i,1}-rawGPS.x{i,1}(1))*1000;
		%rawGPS.y{i,1}=(rawGPS.y{i,1}-rawGPS.y{i,1}(1))*1000;
		rawGPS.z{i,1}=(rawGPS.z{i,1}-rawGPS.z{i,1}(1))*1000;
		%rawGPS.sigx{i,1}=(rawGPS.sigx{i,1})*1000;
		%rawGPS.sigy{i,1}=(rawGPS.sigy{i,1})*1000;
		rawGPS.sigz{i,1}=(rawGPS.sigz{i,1})*1000;

		
	end

		%plotting GPS locations on map
		
		%figure(1);
		%s2=scatter(GPSlong,GPSlat,'^k','MarkerFaceColor','k');
		%legend([s1,s2],'CRMS Stations','GPS Stations')
		%hold off;

		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plotting spatially close CRMS and GPS stations to compare data
%

	% converting lat/longs into UTM coordinates for CRMS and GPS stations
	% this allows us to determine the closest corresponding station using min. distance
	[X,Y]=ll2utm(stationlat,stationlon,15); 
	[X2,Y2]=ll2utm(GPSlat,GPSlong,15);
	[coastx,coasty]=ll2utm(coastlat,coastlon,15);
	[parishx,parishy]=ll2utm(parishlat,parishlon,15);

		X=X/1e3; Y=Y/1e3;
		X2=X2/1e3;	Y2=Y2/1e3;
		coastx=coastx/1e3;	coasty=coasty/1e3;
		parishx=parishx/1e3; parishy=parishy/1e3;         

	% plot check to make sure conversion worked correctly for both station types
	figure(5); clf;
	hold on;
	plot(coastx,coasty,'color',[1 1 1]*0.5); axis equal;
	plot(parishx,parishy,'color',[1 1 1]*0.5);
	plot(X,Y,'.b'); axis equal;
	plot(X2,Y2,'^k','MarkerFaceColor','k');
	hold off;
	xlim([250,950]); ylim([3150,3700]);
	xlabel('Easting (km)')
	ylabel('Northing (km)')
	title('Station Positions in UTM Coordinates')


	stop

	% creating a grid 
	[Ycrms,Ygps]=meshgrid(Y,Y2);
	[Xcrms,Xgps]=meshgrid(X,X2);
	
	% plot check to visualize data
	%figure(6);clf;
	%subplot(221); imagesc(Xgps); colorbar;
	%subplot(222); imagesc(Xcrms); colorbar;
	%subplot(223); imagesc(Ygps); colorbar;
	%subplot(224); imagesc(Ycrms); colorbar;

	% takes the x and y position of each station to calculate the closest corresponding station
	D=sqrt((Xgps-Xcrms).^2+(Ygps-Ycrms).^2);

	% plot check
	%figure(7); clf;
	%imagesc(D); colorbar;

	% uses the minimum distance of each GPS station relative to a CRMS station to create
	% cell arrays with distance and index
	[mindistance,imindistance]=min(D);

	% plots CRMS station with color representing index of closest GPS station
	%figure(8); clf;
	%scatter(X,Y,50,imindistance,'filled'); colorbar;

	%figure(9); clf;
	%scatter(X,Y,50,mindistance,'filled'); colorbar;
	
 	for i=1:numel(stationset)

		figure(10); clf
		set(gcf,'Position',[488.2000   65.8000  635.2000  696.0000])
		subplot(4,1,[1,2]);
		plot(X,Y,'.');
		hold on;
		plot(coastx,coasty,'color',[1 1 1]*0.5); axis equal;
		plot(parishx,parishy,'color',[1 1 1]*0.5);
		plot(X2,Y2,'^k','MarkerFaceColor','k','MarkerSize',2);
		plot(X2(imindistance(i)),Y2(imindistance(i)),'^g','MarkerFaceColor','g');
		plot(X(i),Y(i),'rs','MarkerFaceColor','r');
		xlim([250,950]); ylim([3150,3700]);
		xlabel('Easting (km)')
		ylabel('Northing (km)')
		title([stationset{i} ' vs. GPS ' GPSnames{imindistance(i)} ' ' num2str(mindistance(i),'%0.1f') ' km apart']);
		
		subplot(4,1,[3,4]);
		b=plot(rawGPS.t{imindistance(i)},rawGPS.z{imindistance(i)},'.b'); datetick;
		hold on;
		a=plot(T{i},Hrel{i},'-ro'); datetick
		grid on;
		hold off;
		xlabel('Time')
		ylabel('Vertical Displacement (mm)')
		legend([a,b],'CRMS Data', 'GPS Data','Location','northwest')

		saveas(gcf,['CRMS_Figures/' stationset{i}],'bmp');
		%saveas(gcf,['CRMS_Figures/' stationset{i}],'epsc2');

	end


