#!/usr/local/bin/octave --no-window-system
%#!/usr/bin/octave -qf
%GPSTEC_SH   GPS-TEC data processor (with shebang)
% GPSTEC_SH processes the ground-based GPS observables in RINEX 1.00 or later
% and outputs the vertical TEC and the corresponding longitudes and latitudes
% at the ionospheric points for the specified site and period.
%
% SYNTAX:
%
%   gpstec_sh.m SITE Y D H R MINEL
%
% estimates TEC for site SITE in 4-char, year Y, day-of-year D and hour H in R-
% second sampling above the minimum (cutoff) elevation angle of satellites,
% MINEL, in degrees.
%
%   gpstec_sh.m SITE Y D H1 H2 R MINEL
%
% specifies the period in hours [H1 H2].
%
% For example,
%
%   gpstec_sh.m lpgs 2014 306 0 1 10
%
% estimates TEC for site 'lpgs', date 2014.306, hour 0 (i.e., 000000-005959 GPS
% time), 1-s sampling, above 10-degree elevation angle.
%
%   gpstec_sh.m tash 2014 306 0 24 30 20
%
% estimates TEC for site 'tash', date 2014.306, hours 0 to 24 (i.e., 000000-
% 235930 GPS time), 30-s sampling, above 20-degree elevation angle.
%
% Notice:  Non-GPS data will be discarded.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/21
%   VERSION	:	7.2 2016/12/15

% Set default options to be MATLAB-compatible
if exist('default_save_options','builtin')
	default_save_options('-7'); % Set to MATLAB(R) v7 binary data format
end
if exist('save_default_options','builtin')
	save_default_options('-7');
end

% Handle the input arguments
arg  = argv();
Site = arg{1};
yr   = str2double(arg{2});
doy  = str2double(arg{3});
if size(arg,1)==6
	hrng = str2double(arg{4});
elseif size(arg,1)==7
	hrng = [str2double(arg{4}) str2double(arg{5})];
else
	error('Too many arguments.');
end
rate   = str2double(arg{end-1});
minel  = str2double(arg{end});

% Parameters
clight = 2.99792458d8; % Light speed in m/s
f1 = 1.57542d9; % L1 frequency in Hz
f2 = 1.2276d9; % L2 frequency in Hz
f12 = f1^2; f22 = f2^2;
ionfac = 40.3082d16; % Ionospheric factor in m/TECU/s^2
m2tecu = f12*f22/(f12-f22)/ionfac; % Conversion factor from meters to TECU
ns2tecu = clight*1e-9*m2tecu; % Conversion factor from nanoseconds to TECU
re = 6371; % Mean radius of earth in km
h1 = 650; % Maximum height of ionosphere in km
h2 = 250; % Minimum height of ionosphere in km
ih = (h1+h2)/2; % Mean ionospheric height in km
nepoch = 3600/rate; % Number epochs in an hour
nsat = 32; % Number of satellites
date_as = datenum(1994,1,31); % Date number of anti-spoofing (A/S) switch on
minobs = 3000; % Minimum number of observables for estimating receiver's DCB
GNSS = 'G'; % GNSS type
ptmin = [4 4.5]; % Period of TEC minimum in hours
km2m = 1e3; % Factor between km and meter
rad = 180. /pi;

% Deal with site and time information
Site_lc = lower(Site); % Lowercase site name
Site_uc = upper(Site); % Uppercase site name
Y = num2str(yr,'%04u'); Y2 = Y(3:4);
D = num2str(doy,'%03u'); % Day of year
R = num2str(rate,'%02u'); % Data sampling in seconds
[~,mon,day] = doy2ymd(yr,doy);
M = num2str(mon,'%02u');
DAY = num2str(day,'%02u');
[w,wd] = ymd2gps(yr,mon,day);
W = num2str(w,'%04u');
WD = int2str(wd);
if size(hrng,2)==1
	hrng = [hrng hrng+1];
end
HR = num2str(hrng(1),'%02u');
DM = num2str((hrng(2)-hrng(1))*60-0.01,'%7.2f'); % Time difference in char-minutes
SD = [Site D];

% Load the interpolated orbit file
OUTFILE1 = ['COD' W WD 'i' R '.mat']; % CODwwwwdiRR.mat
if ~exist(OUTFILE1,'file')
	INFILE2 = OUTFILE1([1:8 12:end]); % CODwwwwd.mat
	if exist(INFILE2,'file')
		sp3cubic(INFILE2,rate); % Create CODwwwwdiRR.mat
	elseif exist([INFILE2(1:8) '.EPH'],'file') % CODwwwwd.EPH
		sp32mat(w,wd); % Create CODwwwwd.mat
		sp3cubic(INFILE2,rate); % Create CODwwwwdiRR.mat
	else
		for i=[115 114 117] % ASCII numbers of 's', 'r' and 'u'
			OUTFILE1 = ['ig' char(i) W WD 'i' R '.mat']; % ig[sru]wwwwdiRR.mat
			if exist(OUTFILE1,'file')
				break;
			else
				INFILE2 = OUTFILE1([1:8 12:end]); % ig[sru]wwwwd.mat
				if exist(INFILE2,'file')
					sp3cubic(INFILE2,rate); % Create ig[sru]wwwwdiRR.mat
					break;
				elseif exist([INFILE2(1:8) '.sp3'],'file') % ig[sru]wwwwd.sp3
					sp32mat(w,wd); % Create ig[sru]wwwwd.mat
					sp3cubic(INFILE2,rate); % Create ig[sru]wwwwdiRR.mat
					break;
				end
			end
		end
	end
end
if ~exist(OUTFILE1,'file')
	error('Cannot find any orbit file.');
end
S = load(OUTFILE1,'xi','yi','zi','years','months','days','hours');

% Trim the orbit data according to the specified period
l = S.years==yr & S.months==mon & S.days==day & S.hours>=hrng(1) & ...
	S.hours<hrng(2);
xi  = S.xi(l,1:min(nsat,size(S.xi,2))); % in km NOTICE: GPS only
yi  = S.yi(l,1:min(nsat,size(S.yi,2)));
zi  = S.zi(l,1:min(nsat,size(S.zi,2)));
clear S;

% Compensate columns if not sufficient
[nr,nc] = size(xi);
if nc<nsat
	xi = [xi NaN(nr,nsat-nc)];
end
[nr,nc] = size(yi);
if nc<nsat
	yi = [yi NaN(nr,nsat-nc)];
end
[nr,nc] = size(zi);
if nc<nsat
	zi = [zi NaN(nr,nsat-nc)];
end

% Deal with the O-file(s)
if rate==30
	OUTFILE3 = [SD '0.mat']; % SSSSDDD0.mat
	if exist(OUTFILE3,'file')
		dflag = true; % Flag (1: daily; 0: other)
		if checkvar(OUTFILE3,'rinex_version')
			load(OUTFILE3,'rinex_version');
			rv = floor(rinex_version);
		end
	else
		INFILE3 = [OUTFILE3(1:9) Y2 'o']; % SSSSDDD0.YYo
		if exist(INFILE3,'file')
			rv = floor(str2double(textread(INFILE3,'%s',1)));
			switch rv
			case {1,2}
				system(['teqc -phc -st ' Y2 M DAY HR '0000 +dm ' DM ' -O.obs L1L2C1P1P2 -O.dec ' R 's ' INFILE3 ' > ' INFILE3 '.t1']); % SSSSDDD0.YYo to SSSSDDD0.YYo.t1
				sflag = obs2mat([INFILE3 '.t1']); % SSSSDDD0.YYo.t1 to SSSSDDD0.mat
				if sflag
					delete([INFILE3 '.t1']);
					warning('The site is in the skip list or no any observation data found.');
					return;
				end
				dflag = true;
				delete([INFILE3 '.t1']); % Remove the temporary file
			case 3
				sflag = obs2mat(INFILE3); % SSSSDDD0.YYo to SSSSDDD0.mat
				if sflag
					warning('The site is in the skip list or no observation data found.');
					return;
				end
				dflag = true;
			otherwise
				error('The RINEX version of the input file should be 2.00 or later.');
			end
		else
			error(['Cannot find ' INFILE3 '.']);
		end
	end
elseif rate==1
	for i=hrng(1):hrng(2)-1
		OUTFILE3 = [SD char(97+i) '.mat']; % SSSSDDDA.mat
		if exist(OUTFILE3,'file')
			dflag = false; % Flag (1: daily; 0: other)
			if checkvar(OUTFILE3,'rinex_version')
				load(OUTFILE3,'rinex_version');
				rv = floor(rinex_version);
			end
		else
			INFILE3 = [OUTFILE3(1:9) Y2 'o']; % SSSSDDDA.YYo
			if exist(INFILE3,'file')
				rv = floor(str2double(textread(INFILE3,'%s',1)));
				switch rv
				case {1,2}
					system(['teqc -phc -O.obs L1L2C1P1P2 -O.dec ' R 's ' INFILE3 ' > ' INFILE3 '.t1']); % SSSSDDDA.YYo to SSSSDDA.YYo.t1
					sflag = obs2mat([INFILE3 '.t1']); % SSSSDDDA.YYo.t1 to SSSSDDDA.mat
					if sflag
						delete([INFILE3 '.t1']);
						warning('The site is in the skip list or no observation data found.');
						return;
					end
					dflag = false;
					delete([INFILE3 '.t1']); % Remove the temporary file
				case 3
					sflag = obs2mat(INFILE3); % SSSSDDD0.YYo to SSSSDDD0.mat
					if sflag
						warning('The site is in the skip list or no any observation data found.');
						return;
					end
					dflag = false;
				otherwise
					error('The RINEX version of the input file should be 2.00 or later.');
				end
			else
				error(['Cannot find any O-files for site ' Site '.']);
			end
		end
	end
else
	error('Cannot handle data with sampling other than 1 or 30 seconds.');
end

% Write site's position to file and calculate the elevation angle
load(OUTFILE3,'Domes','x','y','z','TimeSystem');
OUTFILE4 = [Site_uc '_' Domes '.mat']; % SSSS nnnnn[MS]nnnXX.mat
OUTFILE4(OUTFILE4==' ') = '_'; % SSSS_nnnnn[MS]nnnXX.mat
if exist('x','var') && exist('y','var') && exist('z','var') && ~exist(OUTFILE4,'file')
	save(OUTFILE4,'x','y','z');
	disp(['FILE CREATION: ' OUTFILE4]);
end
el = elevation(Site,xi,yi,zi,minel,Domes);

% Calculate longitudes and latitudes of the ionospheric points
[lat,lon] = ionpt(Site,xi,yi,zi,el,ih,Domes);

% Create anti-spoofing flag (True for A/S on)
if datenum(yr,mon,day)>=date_as
	asflag = true;
else
	asflag = false;
end

% Deal with satellite's DCBs
OUTFILE6 = ['P1P2' Y2 M '_ALL.mat']; % P1P2YYMM_ALL.mat
if ~exist(OUTFILE6,'file')
	INFILE6 = ['DCB/' OUTFILE6(1:12) '.DCB']; % DCB/P1P2YYMM_ALL.DCB
	if exist(INFILE6,'file')
		dcb2mat(INFILE6);
	else
		warning([INFILE6 ' not found.']);
	end
end
gpsdcb = zeros(1,nsat); % Allocate space for GPS DCB
if exist(OUTFILE6,'file')
	load(OUTFILE6,'PRN','dcb');
	if exist('PRN','var') && exist('dcb','var')
		for i=1:size(PRN,1)
			prn = str2double(PRN(i,2:3));
			if PRN(i,1)==GNSS && prn<=nsat
				gpsdcb(prn) = dcb(i);
			end
		end
	end
else
	warning([OUTFILE6 ' not found.']);
end
OUTFILE7 = ['P1C1' Y2 M '_RINEX.mat']; % P1C1YYMM_RINEX.mat
if asflag % A/S on
	if ~exist(OUTFILE7,'file')
		INFILE7 = ['DCB/' OUTFILE7(1:14) '.DCB']; % DCB/P1C1YYMM_RINEX.DCB
		if exist(INFILE7,'file')
			dcb2mat(INFILE7);
		else
			warning([INFILE7 ' not found.']);
		end
	end
	if any(gpsdcb)
		if exist(OUTFILE7,'file')
			load(OUTFILE7,'PRN','dcb');
			if exist('PRN','var') && exist('dcb','var')
				for i=1:size(PRN,1)
					prn = str2double(PRN(i,2:3));
					if PRN(i,1)==GNSS && prn<=nsat
						gpsdcb(prn) = gpsdcb(prn)-dcb(i);
					end
				end
			end
		else
			warning([OUTFILE7 ' not found.']);
		end
	end
end

% Read in receiver's DCB
ReceiverName = [Site_uc ' ' Domes]; % Receiver's full name
if exist(OUTFILE6,'file')
	load(OUTFILE6,'PRN','StationName','dcb');
	for i=1:size(PRN,1)
		if PRN(i,1)==GNSS && strcmp(StationName(i,:),ReceiverName),
			recdcb = dcb(i);
		end
	end
end
if asflag && exist('recdcb','var') && exist(OUTFILE7,'file') % Search for P1C1 only if P1P2 is available
	load(OUTFILE7,'PRN','StationName','dcb');
	for i=1:size(PRN,1)
		if PRN(i,1)==GNSS && strcmp(StationName(i,1:4),ReceiverName(1:4)), % NOTICE: The DOMES numbers are missing in P1C1 DCB files and not possible to be compared
			recdcb = recdcb-dcb(i);
		end
	end
end
if ~exist('recdcb','var')
	recdcb = 0;
end

% Solve ambiguity and slant TEC
[ep,stecl,amb] = ambiguity(Site,yr,doy,hrng,el,recdcb,gpsdcb,dflag,rate);
stec = stecl+amb;

% Solve receiver's and/or satellites's DCBs for daily data with sufficient
% observables, or read in STATION.DCB, and then re-estimate the slant TEC if
% applicable
if ~any(gpsdcb) && dflag && nnz(stec)>=minobs
	[recdcb, gpsdcb] = getalldcb(Site_lc,stec,el,ep,xi,yi,zi,Domes);
	l = logical(stec);
	stec = stec+(recdcb+repmat(gpsdcb,size(stec,1),1))*ns2tecu;
	stec(~l) = 0;
end
if ~recdcb
	if dflag && nnz(stec)>=minobs
		recdcb = getrecdcb(Site_lc,stec,el,ep,xi,yi,zi,Domes);
	else
		INFILE8 = 'DCB/STATION.DCB';
		OUTFILE8 = 'STATION.mat';
		if ~exist(OUTFILE8,'file')
			dcb2mat(INFILE8);
		end
		load(OUTFILE8,'PRN','StationName','dcb');
		for i=1:size(PRN,1)
			if PRN(i,1)==GNSS && strcmp(StationName(i,:),ReceiverName),
				recdcb = dcb(i);
			end
		end
	end
	if recdcb
		l = logical(stec);
		stec(l) = stec(l)+recdcb*ns2tecu;
	else
		warning('RDCB not found.  Use 0 instead.');
	end
end

% Get the longitude of the specified site
[rx,ry,rz] = readcrd(Site,Domes); % XYZ in km
[~,rlon] = ecef2lla(rx*km2m,ry*km2m,0);
rlon = rlon*rad;

% Shift all TEC until min(TEC) = 5 if min(TEC) < 5 TECU
ep_lt = ep+rlon/15; % Convert GPS time (almost UT) into local time in hours
l = (ep_lt>ptmin(1) & ep_lt<ptmin(2)) | (ep_lt>ptmin(1)-24 & ...
	ep_lt<ptmin(2)-24) | (ep_lt>ptmin(1)+24 & ep_lt<ptmin(2)+24);
s4 = stec(l,:);
e4 = el(l,:);
l4 = logical(s4) & logical(e4);
S = slant(re,e4(l4),h1,h2);
s4nz = s4(l4); % Non-zero slant TEC
v4nz = s4nz./S; % Non-zero vertical TEC
if ~isempty(v4nz)
	mvtec = median(v4nz); % Median vertical TEC
	if mvtec<5
		offset = -median(s4nz)+5*mean(S); % Slant TEC offset in TECU
		l = logical(stec);
		stec(l) = stec(l)+offset; % Shift TEC until min(TEC) = 5 TECU
		recdcb = recdcb+offset/ns2tecu; % Compensate receiver DCB
	end
end

% Estimate VTEC
[vt,mintec,rlat,rlon] = vtec(Site_lc,doy,ep,stec,el,rate,Domes);

% Synchronize the zero elements in VTEC with other arrays
l = logical(vt);
lon(~l) = 0;
lat(~l) = 0;
el(~l) = 0;

% Write to MAT-file
eval(['v' SD ' = sparse(vt);']);
eval(['s' SD ' = sparse(stec);']);
eval(['o' SD ' = sparse(lon);']);
eval(['a' SD ' = sparse(lat);']);
eval(['e' SD ' = sparse(el);']);
OUTFILE = [SD '.mat']; % SSSSDDD.mat
save(OUTFILE,['v' SD],['s' SD],['o' SD],['a' SD],['e' SD],'ep','recdcb',...
	'gpsdcb','yr','doy','Domes','rlat','rlon','mintec','TimeSystem');
disp(['FILE CREATION: ' OUTFILE]);

% Output ASCII data
OUTFILE = [SD '.dat']; % SSSSDDD.dat
tout = datevec(ep/24);
fid = fopen(OUTFILE,'w');
for i=1:size(vt,1)
	for j=1:nsat
		if vt(i,j)>0
			fprintf(fid,'%2d ',tout(i,4)');
			fprintf(fid,'%2d ',tout(i,5)');
			fprintf(fid,'%2d ',tout(i,6)');
			fprintf(fid,'%s ',Site_lc);
			fprintf(fid,'%2d ',j);
			fprintf(fid,'%8.2f ',lon(i,j));
			fprintf(fid,'%7.2f ',lat(i,j));
			fprintf(fid,'%7.2f \n',vt(i,j));
		end
	end
end
fclose(fid);
disp(['FILE CREATION: ' OUTFILE]);

% Write to the log file
LOGFILE = [SD '.log'];
if exist(LOGFILE,'file')
        delete(LOGFILE); % Remove log file
end
fid = fopen(LOGFILE,'w');
fprintf(fid,'GPS site name: %s\n',[Site ' ' Domes]);
fprintf(fid,'Julian date: %s\n',D);
fprintf(fid,'Site location: %6.2f N, %7.2f E\n',rlat,rlon);
fprintf(fid,'DCB of the receiver: %g ns\n',recdcb);
fprintf(fid,'Minimum of VTEC: %6.2f TECU\n',mintec);
fclose(fid);
disp(['FILE CREATION: ' LOGFILE]);

