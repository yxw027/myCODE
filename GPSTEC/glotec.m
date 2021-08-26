function glotec(Site,yr,doy,hrng,rate,minel,ippalt,inpath,outpath)
%GLOTEC   GLONASS-TEC data processor
% GLOTEC('SITE',Y,D,H,R,MINEL) processes the ground-based GLONASS observables
% in RINEX 1.00 or later and outputs the vertical TEC and the corresponding
% longitudes and latitudes at the ionospheric points for site SITE in 4-char,
% year Y, day-of-year D and hour H in R-second sampling above the minimum
% (cutoff) elevation angle of satellites, MINEL, in degrees.
%
% GLOTEC('SITE',Y,D,[H1 H2],R,MINEL) specifies the period in hours [H1 H2].
%
% For example,
%
%   GLOTEC('lpgs',2014,306,0,1,10)
%
% estimates TEC for site 'lpgs', date 2014.306, hour 0 (i.e., 000000-005959 for
% the specified time system), 1-s sampling, above 10-degree elevation angle.
%
%   GLOTEC('tash',2014,306,[0 24],30,20)
%
% estimates TEC for site 'tash', date 2014.306, hours 0 to 24 (i.e., 000000-
% 235930 for the specified time system), 30-s sampling, above 20-degree
% elevation angle.
%
% Notice:  Non-GLONASS data will be discarded.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/12/03
%   VERSION	:	5.0 2017/02/24

error(nargchk(9,9,nargin));

% Set default options to be MATLAB-compatible
if exist('default_save_options','builtin')
	default_save_options('-7'); % Set to MATLAB(R) v7 binary data format
end
if exist('save_default_options','builtin')
	save_default_options('-7');
end

% Parameters
clight = 2.99792458d8; % Light speed in m/s
f1o = 1.602d9; % Base frequency of L1 in Hz
f1i = 0.5625d6; % Frequency increment of L1 in Hz
f2o = 1.246d9; % Base frequency of L2 in Hz
f2i = 0.4375d6; % Frequency increment of L2 in Hz
ionfac = 40.3082d16; % Ionospheric factor in m/TECU/s^2
re = 6371; % Mean radius of earth in km
h1 = 650; % Maximum height of ionosphere in km
h2 = 250; % Minimum height of ionosphere in km
ih = (h1+h2)/2; % Mean ionospheric height in km
nepoch = 3600/rate; % Number epochs in an hour
nsat = 26; % Number of satellites
minobs = 3000; % Minimum number of observables for estimating receiver's DCB
GNSS = 'R'; % GNSS type
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


tmpdir =['TMP/glonass/' Y '.' D '/'];
obsmatdir =['OBS/glonass/' Y '.' D '/'];
ORBPATH = ['ORB/'];
if ~exist(outpath)
   mkdir(outpath);
end
if ~exist(tmpdir)
   mkdir(tmpdir);
end
if ~exist(obsmatdir)
   mkdir(obsmatdir);
end
if ~exist(ORBPATH)
   mkdir(ORBPATH);
end

% Load the interpolated orbit file
OUTFILE1 = ['COD' W WD 'i' R '.mat']; % CODwwwwdiRR.mat
if ~exist([ORBPATH OUTFILE1],'file')
	INFILE2 = OUTFILE1([1:8 12:end]); % CODwwwwd.mat
	if exist([ORBPATH INFILE2],'file')
		sp3cubic([ORBPATH INFILE2],rate); % Create CODwwwwdiRR.mat
	elseif exist([ORBPATH INFILE2(1:8) '.EPH'],'file') % CODwwwwd.EPH
		sp32mat(w,wd); % Create CODwwwwd.mat
		sp3cubic([ORBPATH INFILE2],rate); % Create CODwwwwdiRR.mat
	else
		error('Cannot find any orbit file.');
	end
end
S = load([ORBPATH OUTFILE1],'xi','yi','zi','years','months','days','hours');

% Trim the orbit data according to the specified period
l = S.years==yr & S.months==mon & S.days==day & S.hours>=hrng(1) & ...
	S.hours<hrng(2);
xi  = S.xi(l,33:end); % in km NOTICE: GLONASS only
yi  = S.yi(l,33:end);
zi  = S.zi(l,33:end);
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
	if exist([obsmatdir OUTFILE3],'file')
		dflag = true; % Flag (1: daily; 0: other)
		if checkvar([obsmatdir OUTFILE3],'rinex_version')
			load([obsmatdir OUTFILE3],'rinex_version');
			rv = floor(rinex_version);
		end
	else
		INFILE3 = [OUTFILE3(1:9) Y2 'o']; % SSSSDDD0.YYo
		if exist([inpath INFILE3],'file')
			rv = floor(str2double(textread([inpath INFILE3],'%s',1)));
			switch rv
			case {1,2}
                tmpfile = [ OUTFILE3(1:9) Y2 'o.t1'];
				system(['teqc -phc -st ' Y2 M DAY HR '0000 +dm ' DM ' -O.obs L1L2C1P2 -O.dec ' R 's ' inpath INFILE3 ' > ' pwd '/' tmpfile]); % SSSSDDD0.YYo to SSSSDDD0.YYo.t1  NOTICE: P1 is discarded
				sflag = obs2mat_m(tmpfile,GNSS); % SSSSDDD0.YYo.t1 to SSSSDDD0.mat
				if sflag
					delete(tmpfile);
					warning('The site is in the skip list or no any observation data found.');
					return;
				end
				dflag = true;
                delete(tmpfile); % Remove the temporary file
                movefile(OUTFILE3,obsmatdir);
			case 3
                copyfile([inpath INFILE3],pwd); 
				sflag = obs2mat_m(INFILE3,GNSS); % SSSSDDD0.YYo to SSSSDDD0.mat
				if sflag
					warning('The site is in the skip list or no observation data found.');
					return;
				end
				flag = true;
                movefile(OUTFILE3,obsmatdir);
                delete(INFILE3);
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
		if exist([obsmatdir OUTFILE3],'file')
			dflag = false; % Flag (1: daily; 0: other)
			if checkvar([obsmatdir OUTFILE3],'rinex_version')
				load([obsmatdir OUTFILE3],'rinex_version');
				rv = floor(rinex_version);
			end
		else
			INFILE3 = [OUTFILE3(1:9) Y2 'o']; % SSSSDDDA.YYo
			if exist([inpath INFILE3],'file')
				rv = floor(str2double(textread([inpath INFILE3],'%s',1)));
				switch rv
				case {1,2}
                    tmpfile = [ OUTFILE3(1:9) Y2 'o.t1'];
					system(['teqc -phc -O.obs L1L2C1P2 -O.dec ' R 's ' inpath INFILE3 ' > ' pwd '/' tmpfile]); % SSSSDDDA.YYo to SSSSDDA.YYo.t1  NOTICE: P1 is discarded
					sflag = obs2mat_m(tmpfile,GNSS); % SSSSDDDA.YYo.t1 to SSSSDDDA.mat
					if sflag
						delete(tmpfile);
						warning('The site is in the skip list or no observation data found.');
						return;
					end
					dflag = false;
					delete(tmpfile); % Remove the temporary file
                    movefile(OUTFILE3,obsmatdir);
				case 3
                    copyfile([inpath INFILE3],pwd);
					sflag = obs2mat_m(INFILE3,GNSS); % SSSSDDD0.YYo to SSSSDDD0.mat
					if sflag
						warning('The site is in the skip list or no any observation data found.');
						return;
					end
					dflag = false;
                    movefile(OUTFILE3,obsmatdir);
                    delete(INFILE3);
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
load([obsmatdir OUTFILE3],'Domes','x','y','z','TimeSystem');
OUTFILE4 = [Site_uc '_' Domes '.mat']; % SSSS nnnnn[MS]nnnXX.mat
OUTFILE4(OUTFILE4==' ') = '_'; % SSSS_nnnnn[MS]nnnXX.mat
if exist('x','var') && exist('y','var') && exist('z','var') && ~exist(OUTFILE4,'file')
	save(OUTFILE4,'x','y','z');
	disp(['FILE CREATION: ' OUTFILE4]);
end
el = elevation_r(Site,xi,yi,zi,minel,Domes);

% Calculate longitudes and latitudes of the ionospheric points
[lat,lon] = ionpt_r(Site,xi,yi,zi,el,ippalt,Domes);

% Deal with satellite's DCBs
OUTFILE6 = ['P1P2' Y2 M '_ALL.mat']; % P1P2YYMM_ALL.mat
if ~exist(['DCB/' OUTFILE6],'file')
	INFILE6 = ['DCB/' OUTFILE6(1:12) '.DCB']; % DCB/P1P2YYMM_ALL.DCB
	if exist(INFILE6,'file')
		dcb2mat(INFILE6);
	else
		warning([INFILE6 ' not found.']);
	end
end
glodcb = zeros(1,nsat); % Allocate space for GLO DCB
if exist(['DCB/' OUTFILE6],'file')
	load(['DCB/' OUTFILE6],'PRN','dcb');
	if exist('PRN','var') && exist('dcb','var')
		for i=1:size(PRN,1)
			prn = str2double(PRN(i,2:3));
			if PRN(i,1)==GNSS && prn<=nsat
				glodcb(prn) = dcb(i);
			end
		end
	end
else
	warning([OUTFILE6 ' not found.']);
end
OUTFILE7 = ['P1C1' Y2 M '_RINEX.mat']; % P1C1YYMM_RINEX.mat
if ~exist(['DCB/' OUTFILE7],'file')
	INFILE7 = ['DCB/' OUTFILE7(1:14) '.DCB']; % DCB/P1C1YYMM_RINEX.DCB
	if exist(INFILE7,'file')
		dcb2mat(INFILE7);
	else
		warning([INFILE7 ' not found.']);
	end
end
if any(glodcb)
	if exist(['DCB/' OUTFILE7],'file')
		load(['DCB/' OUTFILE7],'PRN','dcb');
		if exist('PRN','var') && exist('dcb','var')
			for i=1:size(PRN,1)
				prn = str2double(PRN(i,2:3));
				if PRN(i,1)==GNSS && prn<=nsat
					glodcb(prn) = glodcb(prn)-dcb(i);
				end
			end
		end
	else
		warning([OUTFILE7 ' not found.']);
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
if exist('recdcb','var') && exist(OUTFILE7,'file') % Search for P1C1 only if P1P2 is available
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

% Fetch frequency table and create conversion factor table
switch rv
case 2
	nav2mat(['brdc' D '0.' Y2 'g']); % brdcDDD0.YYg
	INFILE = ['brdc' D 'g.mat'];
	load(['EPM/' INFILE],'sid','freqn');
case 3
	load(OUTFILE3,'sid','freqn');
end
k = zeros(1,nsat);
for i=1:nsat
	f = find(sid==i);
	if f
		k(i) = freqn(f(1));
	end
end
f1 = f1o+f1i*k; % L1 frequencies in Hz
f2 = f2o+f2i*k; % L2 frequencies in Hz
f12 = f1.^2; f22 = f2.^2;
m2tecu = f12.*f22./(f12-f22)/ionfac; % Conversion factor from meters to TECU
ns2tecu = clight*1e-9*m2tecu; % Conversion factor from nanoseconds to TECU

% Solve ambiguity and slant TEC
[ep,stecl,amb] = ambiguity_r(Site,yr,doy,hrng,el,recdcb,glodcb,dflag,rate);
stec = stecl+amb;

% Solve receiver's DCB for daily data with sufficient observables or read in
% STATION.DCB, and then re-estimate the slant TEC if applicable
if ~recdcb
	if dflag && nnz(stec)>=minobs
		recdcb = getrecdcb_r(Site_lc,yr,doy,stec,el,ep,xi,yi,zi,sid,freqn,Domes);
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
		offset = repmat(recdcb*ns2tecu,size(stec,1),1); % recdcb in TECU
		stec(l) = stec(l)+offset(l);
	else
		warning('RDCB not found.  Use 0 instead.');
	end
end

% Get the longitude of the specified site
[rx,ry,rz] = readcrd(Site,Domes); % XYZ in km
[~,rlon] = ecef2lla(rx*km2m,ry*km2m,0);
rlon = rlon*rad;

% Shift all TEC until min(TEC) = 5 if min(TEC) < 5 TECU
ep_lt = ep+rlon/15; % Convert the given time (almost UT) into local time in hours
l = (ep_lt>ptmin(1) & ep_lt<ptmin(2)) | (ep_lt>ptmin(1)-24 & ...
	ep_lt<ptmin(2)-24) | (ep_lt>ptmin(1)+24 & ep_lt<ptmin(2)+24);
s4 = stec(l,:);
e4 = el(l,:);
l4 = logical(s4) & logical(e4);
S = slant(re,e4(l4),h1,h2);
s4nz = s4(l4); % Non-zero slant TEC
v4nz = s4nz./S; % Non-zero vertical TEC
if any(v4nz)
	mvtec = median(v4nz); % Median vertical TEC
	if mvtec<5
		offset = -median(s4nz)+5*mean(S); % Slant TEC offset in TECU
		l = logical(stec);
		stec(l) = stec(l)+offset; % Shift TEC until min(TEC) = 5 TECU
		recdcb = recdcb+offset/median(ns2tecu); % Compensate receiver DCB
	end
end

% Estimate VTEC
[vt,mintec,rlat,rlon] = vtec(Site_lc,doy,ep,stec,el,rate,ippalt,Domes);

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
save([outpath OUTFILE],['v' SD],['s' SD],['o' SD],['a' SD],['e' SD],'ep','recdcb',...
	'glodcb','yr','doy','Domes','rlat','rlon','mintec','TimeSystem');
disp(['FILE CREATION: ' OUTFILE]);

% Write to the log file
LOGFILE = [tmpdir SD '.log'];
if exist(LOGFILE,'file')
	delete(LOGFILE); % Remove log file
end
fid = fopen(LOGFILE,'w');
fprintf(fid,'GLONASS site name: %s\n',[Site ' ' Domes]);
fprintf(fid,'Julian date: %s\n',D);
fprintf(fid,'Site location: %6.2f N, %7.2f E\n',rlat,rlon);
fprintf(fid,'DCB of the receiver: %g ns\n',recdcb);
fprintf(fid,'Minimum of VTEC: %6.2f TECU\n',mintec);
fclose(fid);
disp(['FILE CREATION: ' LOGFILE]);
movefile(OUTFILE4,tmpdir);

