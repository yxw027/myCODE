function sp32mat(w,wd)
% SP32MAT  Convert an SP3 file into a MAT-file.
% SP32MAT(W,WD) converts an SP3 file, named ig[rsu]wwwwd.sp3 or CODwwwwd.EPH
% into a MAT-file for GPS week W, and day of week WD.  Versions a and c are
% both acceptable while GLONASS slot numbers will be shifted by 32 if any,
% i.e., data for R01 are shifted to column 33, R02 to column 34, etc.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/12
%   VERSION	:	3.0 2017/03/08

narginchk(2,2);

% Parameters
W = num2str(w,'%04u');
WD = num2str(wd,'%01u');
ORBPATH = ['ORB/'];
% Load the SP3 file
INFILE = ['COD' W WD '.EPH']; % CODWWWWD.EPH
fid = fopen([ORBPATH INFILE],'r');
if fid==-1
	INFILE = ['igs' W WD '.sp3']; % igsWWWWD.sp3
	fid = fopen([ORBPATH INFILE],'r');
	if fid==-1
		INFILE(3) = 'r'; % igrWWWWD.sp3
		fid = fopen([ORBPATH INFILE],'r');
		if fid==-1
			INFILE(3) = 'u'; % iguWWWWD.sp3
			fid = fopen([ORBPATH INFILE],'r');
			if fid==-1
				error(['Cannot find any SP3 file for GPS week ' W ', day of the week ' WD '.']);
			end
		end
	end	
end

% Read in the first line of the header
Version    = fscanf(fid,'%2c',1); % Version (#a or #c)
PV         = fscanf(fid,'%c',1); % Position and/or velocity
year       = int16(fscanf(fid,'%4d',1)); % Year start
month      = int8(fscanf(fid,'%2d',1)); % Month start
day        = int8(fscanf(fid,'%2d',1)); % Day of month start
hour       = int8(fscanf(fid,'%2d',1)); % Hour start
minute     = int8(fscanf(fid,'%2d',1)); % Minute start
second     = single(fscanf(fid,'%11f',1)); % Second start
nepoch     = int32(fscanf(fid,'%7d',1)); % Number of epochs
Data_used  = fscanf(fid,'%5c',1); % Data used
fscanf(fid,'%c',1); % Skip 1 blank
Coordinate = fscanf(fid,'%5c',1); % Coordinate system
fscanf(fid,'%c',1); % Skip 1 blank
Orbit_type = fscanf(fid,'%3c',1); % Orbit type
fscanf(fid,'%c',1); % Skip 1 blank
Agency     = fscanf(fid,'%4c',1); % Agency
fgetl(fid); % Skip the rest part of the line
 
% Read in the second line of the header
fscanf(fid,'%2c%c',2); % skip symbols
gpsweek    = int16(fscanf(fid,'%4d',1)); % GPS week
sec_week   = single(fscanf(fid,'%g',1)); % seconds of week
epo_intv   = single(fscanf(fid,'%g',1)); % epoch interval in seconds
mjd        = uint16(fscanf(fid,'%5d',1)); % modified Julian day start (max 65535)
fday       = single(fscanf(fid,'%g',1)); % fractional day
fgetl(fid); % Skip the rest part of the line

% Read in the third line of the header
fscanf(fid,'%2c%2c',2); % Skip symbols and 2 blanks
nsat = fscanf(fid,'%2d',1); % Number of satellites
fscanf(fid,'%3c',1); % Skip 3 blanks
id = fscanf(fid,'%c%2d',min(nsat*2,34)); % Scan 17 or less satellite IDs (the originals can be ' nn', 'Gnn' or 'Rnn')
gnss = id(1:2:end-1); % GNSS type, i.e., 'G' or 'R' in ASCII
sid = id(2:2:end); % Take the satellite IDs only (omit G and R)
fgetl(fid); % Skip the rest part of the line

% Read in the fourth line of the header
if nsat>17
	fscanf(fid,'%2c%7c',2); % Skip symbols and 7 blanks
	id = fscanf(fid,'%c%2d',min((nsat-17)*2,34)); % Scan the second rows of PRNs or slot numbers
	gnss = [gnss; id(1:2:end-1)];
	sid = [sid; id(2:2:end)];
end
fgetl(fid); % Skip all or the rest part of the line

% Read in the fifth line of the header
if nsat>34
	fscanf(fid,'%2c%7c',2);
	id = fscanf(fid,'%c%2d',min((nsat-34)*2,34));
	gnss = [gnss; id(1:2:end-1)];
	sid = [sid; id(2:2:end)];
end
fgetl(fid); % Skip all or the rest part of the line

% Read in the sixth line of the header
if nsat>51
	fscanf(fid,'%2c%7c',2);
	id = fscanf(fid,'%c%2d',min((nsat-51)*2,34));
	gnss = [gnss; id(1:2:end-1)];
	sid = [sid; id(2:2:end)];
end
fgetl(fid); % Skip all or the rest part of the line

% Shift data for GLONASS R01 to Column 33, R02 to Column 34, etc.
l = gnss==82;
sid(l) = sid(l)+32;

% Allocate space for arrays
years = zeros(nepoch,1,'int16'); % Array of year start
months = zeros(nepoch,1,'int8'); % Array of month start
days = months; hours = months; minutes = months; % Arrays of day-of-month start, hour and minute start
seconds = zeros(nepoch,1,'single'); % Array of second start
Pv = char(32*ones(nepoch,32)); % Char array of position/velocity (char(32)=' ')
x = NaN(nepoch,max(sid)); % Array for x-coordinate in km
y = x; z = x; t = x; % Arrays for y-, z-coordinates and clock correction
if PV=='V'
	xd = x; yd = y; zd = z; % Velocity arrays for XYZ in dm/s
	td = t; % Clock rate change
end
 
% Skip the rest part of the header
for i=7:22
	fgetl(fid);
end

% Read in the body of the file
for i=1:nepoch
	fscanf(fid,'%2c%c',2); % Skip symbols
	years(i)   = fscanf(fid,'%4d',1); % Year start
	months(i)  = fscanf(fid,'%2d',1); % Month start
	days(i)    = fscanf(fid,'%2d',1); % Day of month start
	hours(i)   = fscanf(fid,'%2d',1); % Hour start
	minutes(i) = fscanf(fid,'%2d',1); % Minute start
	seconds(i) = fscanf(fid,'%g',1); % Second start
	fgetl(fid); % Skip the rest part of the line
	for j=1:nsat
		fscanf(fid,'%c%3c',2); % Skip symbol (P or V) and vehicle ID
		x(i,sid(j))  = fscanf(fid,'%g',1); % X-coordinate in km
		y(i,sid(j))  = fscanf(fid,'%g',1); % Y-coordinate in km
		z(i,sid(j))  = fscanf(fid,'%g',1); % Z-coordinate in km
		t(i,sid(j))  = fscanf(fid,'%g',1); % Clock correction in microsecs
		fgetl(fid); % Skip the rest part of the line
		if PV=='V'
			fscanf(fid,'%c%3c',2); % Skip symbol (V) and vehicle ID
			xd(i,sid(j)) = fscanf(fid,'%g',1); % X-velocity in dm/s
			yd(i,sid(j)) = fscanf(fid,'%g',1); % Y-velocity in dm/s
			zd(i,sid(j)) = fscanf(fid,'%g',1); % Z-velocity in dm/s
			td(i,sid(j)) = fscanf(fid,'%g',1); % Clock rate-change
			fgetl(fid); % Skip the next line
		end
	end
end
fclose(fid);
 
% Replace the bad data with NaNs
t(t==999999.999999) = NaN;

% Convert data type
%sid = int8(sid);

% Save to file
OUTFILE = [ORBPATH INFILE(1:8) '.mat'];
switch PV
case 'P'
	save(OUTFILE,'x','y','z','t','years','months','days','hours',...
		'minutes','seconds','PV','year','month','day','hour',...
		'minute','second','nepoch','Data_used','Coordinate',...
		'Orbit_type','Agency','gpsweek','sec_week','epo_intv',...
		'mjd','fday','nsat','Version');
case 'V'
	save(OUTFILE,'x','y','z','t','xd','yd','zd','td','years','months',...
		'days','hours','minutes','seconds','PV','year','month',...
		'day','hour','minute','second','nepoch','Data_used',...
		'Coordinate','Orbit_type','Agency','gpsweek','sec_week',...
		'epo_intv','mjd','fday','nsat','Version');
otherwise
	error('Cannot find Pos or Vel Flag.');
end
disp(['FILE CREATION: ' OUTFILE]);

