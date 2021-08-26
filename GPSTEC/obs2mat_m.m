function sflag = obs2mat_m(INFILE,GNSS)
%OBS2MAT_M  Convert OBS file into MAT-file for GLONASS, Galileo and BDS.
% FLAG = OBS2MAT_M('SSSSDDDN.YYo','GNSS') converts an OBS file, SSSSDDDN.YYo,
% in RINEX 2.00 or later into a MAT-file for the specified GNSS type (R:
% GLONASS; E: Galileo; C: BeiDou).  FLAG indicates if the file is in the skip
% list or no any observation data found (1: true; 0: false).  For example,
%
%    FLAG = OBS2MAT_M('tacc070a.11o','R')
%
% converts 'tacc070a.11o' into 'tacc070a.mat' and returns 0 if not skipped and
% observation data available; otherwise, returns 1 (no file output).
%
% Notice:  Non-GLONASS, Galileo, BeiDou data will be discarded.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/12/03
%   VERSION	:	3.0 2016/12/22

error(nargchk(2,2,nargin));

% Parameters
Site_uc = upper(INFILE(1:4)); % Site name in uppercase
doy = str2double(INFILE(5:7)); % Day of year
A = INFILE(8); % Alphanumeric flag in char
yr = str2double(INFILE(10:11)); % 2-digit year
m2km = 1d-3; % Conversion factor between m and km
bnd_pr.r =[1.5e7 2.5e7]; % Boundary for GLO code observables in meters
bnd_pr.e =[1.5e7 3e7]; % Boundary for GAL code observables in meters
bnd_pr.c =[1.5e7 4e7]; % Boundary for BDS code observables in meters
sflag = true; % Preset flag as true (no file output yet)
if isunix || ismac
	CMD = 'grep -v ';
elseif ispc
	CMD = 'findstr /V /C:';
else
	error('For UNIX-like, Mac OS X and Windows systems only.');
end

% Check the RINEX version
if exist(INFILE,'file')
	rinex_version = str2double(textread(INFILE,'%s',1));
	rv = floor(rinex_version);
else
	warning('File not found.');
	return;
end
if rv==2 || rv==3
	if rinex_version==3
		warning('Data may contain phase cycle shifts.');
	end
else
	warning('The RINEX version of the input file should be 2.00 or later.');
	return;
end

% Check the GNSS type
switch GNSS
case 'R', GNSSfn = 'GLONASS'; nsat = 26;
case 'E', GNSSfn = 'Galileo'; nsat = 30;
case 'C', GNSSfn = 'BeiDou'; nsat = 35;
otherwise,
	warning('GNSS type can only be ''R'', ''E'', or ''C''.');
	return;
end

% Pre-processing
if rv==2
	switch GNSS
	case 'R'
		system(['teqc -G -E -S ' INFILE ' > ' INFILE '.t1']); % SSSSDDDA.YYo to SSSSDDDA.YYo.t1
	case 'E'
		system(['teqc -G -R -S ' INFILE ' > ' INFILE '.t1']); % SSSSDDDA.YYo to SSSSDDDA.YYo.t1
	case 'C'
		system(['teqc -G -R -E -S -J +relax ' INFILE ' > ' INFILE '.t1']); % SSSSDDDA.YYo to SSSSDDDA.YYo.t1
	end
%	pause(3); % Make sure teqc processing complete
	system([CMD '"                            4  1" ' INFILE '.t1 | ' CMD '"^ .*\.0000000  .  .$" | ' CMD '"COMMENT" > ' INFILE '.t2']); % SSSSDDDA.YYo.t1 to SSSSDDDA.YYo.t2
%	pause(3); % Make sure grep processing complete
	delete([INFILE '.t1']); % Remove the temporary file
end

% Open file
switch rv
case 2, [fid, Msg] = fopen([INFILE '.t2'],'r'); % Open SSSSDDDA.YYo.t2
case 3, [fid, Msg] = fopen(INFILE,'r'); % Open SSSSDDDA.YYo
end
if fid==-1
	warning(Msg);
	return;
end

% Read in the header of the file
S = fgetl(fid);
switch rv
case 2,
	if size(S,2)>=80 && strcmp(S(61:80),'RINEX VERSION / TYPE')
		if S(21)~='O' % File type
			fclose(fid);
			warning('Not observation file input.');
			delete([INFILE '.t2']);
			return;
		end
	else
		fclose(fid);
		warning('RINEX version / type not found.');
		delete([INFILE '.t2']);
		return;
	end
	while ~feof(fid) && ~strcmp(S(61:73),'END OF HEADER')
		S = fgetl(fid);
		S(size(S,2)+1:80) = char(32); % Compensate the tail with whitespaces if necessary
		switch S(61:79)
		case 'MARKER NUMBER      '
			Marker_number = deblank(S(1:11));
			n = size(Marker_number,2);
			if n>=9 && n<11
				Domes = [Marker_number char(32*ones(1,11-n))];
			else
				Domes = char(32*ones(1,11)); % All whitespaces
			end
		case 'APPROX POSITION XYZ'
			x = str2double(S( 1:14))*m2km; % X-coordinate in km (ECF)
			y = str2double(S(15:28))*m2km; % Y-coordinate in km (ECF)
			z = str2double(S(29:42))*m2km; % Z-coordinate in km (ECF)
		case '# / TYPES OF OBSERV'
			ntype = str2double(S(5:6));
			if ntype
				for i=1:min(ntype,9)
					Type(i,:) = S([5 6]+i*6);
				end
			if ntype>9
				S = fgetl(fid);
				for i=10:min(ntype,18)
					Type(i,:) = S([-49 -48]+i*6);
				end
			if ntype>18
				S = fgetl(fid);
				for i=19:min(ntype,20)
					Type(i,:) = S([-103 -102]+i*6);
				end
			if ntype>20
				fclose(fid);
				warning('Too many observation types to be handle.');
				delete([INFILE '.t2']);
				return;
			end, end, end, else
				fclose(fid);
				warning('No any observation types found.');
				delete([INFILE '.t2']);
				return;
			end
		case 'INTERVAL           '
			intv = str2double(S(1:10));
		case 'TIME OF FIRST OBS  '
			TimeSystem = S(49:51);
		end
	end
case 3
	if size(S,2)>=80 && strcmp(S(61:80),'RINEX VERSION / TYPE')
		if S(21)~='O' % File type
			fclose(fid);
			warning('Not observation file input.');
			return;
		end
		Sat_system = S(41); % Satellite system
		if Sat_system~=GNSS && Sat_system~='M'
			fclose(fid);
			warning(['The file contains no ' GNSSfn ' data.']);
			return;
		end
	else
		fclose(fid);
		warning('RINEX version / type not found.');
		return;
	end
	while ~feof(fid) && ~strcmp(S(61:73),'END OF HEADER')
		S = fgetl(fid);
		S(size(S,2)+1:80) = char(32); % Compensate the tail with whitespaces if necessary
		switch S(61:80)
		case 'MARKER NAME         '
			Marker_name = deblank(S(1:60));
		case 'MARKER NUMBER       '
			Marker_number = deblank(S(1:20));
			n = size(Marker_number,2);
			if n>=9 && n<11
				Domes = [Marker_number char(32*ones(1,11-n))];
			else
				Domes = char(32*ones(1,11));
			end
		case 'APPROX POSITION XYZ '
			x = str2double(S( 1:14))*m2km; % X-coordinate in km (ECF)
			y = str2double(S(15:28))*m2km; % Y-coordinate in km (ECF)
			z = str2double(S(29:42))*m2km; % Z-coordinate in km (ECF)
		case 'SYS / # / OBS TYPES '
			if S(1)==GNSS
				ntype = str2double(S(4:6));
				if ntype
					for i=1:min(ntype,13)
						Type(i,:) = S((4:6)+i*4);
					end
				if ntype>13
					S = fgetl(fid);
					for i=14:min(ntype,26)
						Type(i,:) = S((-48:-46)+i*4);
					end
				if ntype>26
					fclose(fid);
					warning('Too many observation types to be handled.');
					return;
				end, end, else
					fclose(fid);
					warning('No any observation types found.');
					return;
				end
			else
				continue;
			end
		case 'INTERVAL            '
			intv = str2double(S(1:10));
		case 'TIME OF FIRST OBS   '
			TimeSystem = S(49:51);
		case 'GLONASS SLOT / FRQ #'
			nslot = str2double(S(1:3)); % Number of slot numbers
			if nslot
				for i=1:min(nslot,8)
					i7 = i*7;
					sid(i) = str2double(S([-1 0]+i7));
					freqn(i) = str2double(S([2 3]+i7));
				end
			if nslot>8
				S = fgetl(fid);
				for i=9:min(nslot,16)
					i7 = i*7;
					sid(i) = str2double(S([-57 -56]+i7));
					freqn(i) = str2double(S([-54 -53]+i7));
				end
			if nslot>16
				S = fgetl(fid);
				for i=17:min(nslot,24)
					i7 = i*7;
					sid(i) = str2double(S([-113 -112]+i7));
					freqn(i) = str2double(S([-110 -109]+i7));
				end
			if nslot>24
				S = fgetl(fid);
				for i=25:min(nslot,26)
					i7 = i*7;
					sid(i) = str2double(S([-169 -168]+i7));
					freqn(i) = str2double(S([-166 -165]+i7));
				end
			end, end, end, end
		end
	end
end

% Stop processing if no any relevant observables found
if feof(fid)
	fclose(fid);
	warning(['No any relevant observables found in ' INFILE '.']);
	if rv==2, delete([INFILE '.t2']); end
	return;
end

% Create a null DOMES number if not available
if ~exist('Domes','var')
	Domes = char(32*ones(1,11));
end

% Create null XYZ if not available
if ~exist('x','var')
	x = [];
	y = [];
	z = [];
end

% Set default time system as GPS if not specified
if ~exist('TimeSystem','var') || strcmp(TimeSystem,'   ')
	TimeSystem = 'GPS';
end

% Stop processing if the station and date are in the skip list
ReceiverName = [Site_uc ' ' Domes]; % Receiver full name
if skipcheck(ReceiverName,yr,doy) % Skip flag (1: to skip; 0: not to skip)
	fclose(fid);
	warning(['Processing terminated as ' ReceiverName ' is in the skip list.']);
	if rv==2, delete([INFILE '.t2']); end
	return;
end

% Allocate space for variables
if A=='0'
	n = 86400. /intv;
else
	n = 3600. /intv;
end
for i=1:ntype
	eval([Type(i,:) ' = zeros(n,nsat);']);
	tmp{i} = zeros(n,nsat); % Create temporary arrays
end

% Read in the body of the file
i = 0;
switch rv
case 2
	while ~feof(fid)
		i = i+1;

		% Deal with the EPOCH / SAT lines
		eps = fscanf(fid,'%u %u %u %u %u %g %u %u',8);
		year(i)       = eps(1);
		month(i)      = eps(2);
		day(i)        = eps(3);
		hour(i)       = eps(4);
		minute(i)     = eps(5);
		second(i)     = eps(6);
		event_flag(i) = eps(7);
		nprn          = eps(8);
		S = fgetl(fid);

		% Check if the PRNs or slot numbers are sufficient (the real number of PRNs may be less than the number of satellites 'nprn' due to unknown issue) 
		if size(S,2)<min(nprn*3,36)
			fclose(fid);
			warning(['Cannot read PRN or slot numbers for ' num2str(year(i),'%02u') '-' num2str(month(i),'%02u') '-' num2str(day(i),'%02u') ' ' num2str(hour(i),'%02u') ':' num2str(minute(i),'%02u') ':' num2str(second(i),'%02u') '.']);
			delete([INFILE '.t2']);
			return;
		end

		m = mod(size(S,2),3);
		if m==1
			S = ['  ' S];
		elseif m==2
			S = [' ' S];
		end
		for j=1:min(nprn,12)
			k = j*3-2;
			if S(k)==GNSS
				prn(j) = str2double(S(k+[1 2]));
			else
				fclose(fid);
				warning(['Cannot handle non-' GNSSfn ' data (number ' S(k+[0 1 2]) ').']);
				delete([INFILE '.t2']);
				return;
			end
		end

		% Read the rest PRNs/slot numbers (second line if more than 12 PRNs)
		if nprn>12
			S = fgetl(fid);
			for j=13:min(nprn,24)
				k = j*3-6;
				if S(k)==GNSS
					prn(j) = str2double(S(k+[1 2]));
				else
					fclose(fid);
					warning(['Cannot handle non-' GNSSfn ' data (number ' S(k+[0 1 2]) ').']);
					delete([INFILE '.t2']);
					return;
				end
			end
		end

		% Stop if more than 24 PRNs/slot numbers
		if nprn>24
			fclose(fid);
			warning('Too many PRN or slot numbers to be handled.');
			delete([INFILE '.t2']);
			return;
		end

		% Deal with the OBSERVATIONS lines
		for j=1:nprn

			% Read observables (first line)
			S = fgetl(fid);
			for k=1:min(ntype,5)
				l = [-15 -2]+k*16;
				if size(S,2)>=l(2)
					val = str2double(S(l(1):l(2)));
					if isfinite(val)
						tmp{k}(i,prn(j)) = val;
					end
				end
			end

			% Read observables (second line if more than 5 data types)
			if ntype>5
				S = fgetl(fid);
				for k=6:min(ntype,10)
					l = [-95 -82]+k*16;
					if size(S,2)>=l(2)
						val = str2double(S(l(1):l(2)));
						if isfinite(val)
							tmp{k}(i,prn(j)) = val;
						end
					end
				end
			end

			% Read observables (third line if more than 10 data types)
			if ntype>10
				S = fgetl(fid);
				for k=11:min(ntype,15)
					l = [-175 -162]+k*16;
					if size(S,2)>=l(2)
						val = str2double(S(l(1):l(2)));
						if isfinite(val)
							tmp{k}(i,prn(j)) = val;
						end
					end
				end
			end

			% Read observables (fourth line if more than 15 data types)
			if ntype>15
				S = fgetl(fid);
				for k=16:min(ntype,20)
					l = [-255 -242]+k*16;
					if size(S,2)>=l(2)
						val = str2double(S(l(1):l(2)));
						if isfinite(val)
							tmp{k}(i,prn(j)) = val;
						end
					end
				end
			end
		end
	end
	fclose(fid);
	pause(0.5); % To make sure file being completely closed
	delete([INFILE '.t2']);
case 3
	while ~feof(fid)

		% Deal with the EPOCH / SAT lines
		[S, count] = fscanf(fid,'%c %u %u %u %u %u %g %u %u',9);
		if S(1)==62 % >
			if count==9
				i = i+1;
				year(i)       = S(2);
				month(i)      = S(3);
				day(i)        = S(4);
				hour(i)       = S(5);
				minute(i)     = S(6);
				second(i)     = S(7);
				event_flag(i) = S(8);
				nprn          = S(9);
				S = fgetl(fid);
				if size(S,2)
					receiver_clk_offset(i) = str2double(S);
				end
			else
				break;
			end
		else
			fclose(fid);
			warning('Cannot find any epoch record.');
			return;
		end

		for j=1:nprn
			S = fgetl(fid);
			n = size(S,2);
			if n>=3 && S(1)==GNSS
				prn = str2double(S(2:3)); % PRN
			if n>=17
				val = str2double(S(4:17));
				if isfinite(val), tmp{1}(i,prn) = val; end
			if n>=33
				val = str2double(S(20:33));
				if isfinite(val), tmp{2}(i,prn) = val; end
			if n>=49
				val = str2double(S(36:49));
				if isfinite(val), tmp{3}(i,prn) = val; end
			if n>=65
				val = str2double(S(52:65));
				if isfinite(val), tmp{4}(i,prn) = val; end
			if n>=81
				val = str2double(S(68:81));
				if isfinite(val), tmp{5}(i,prn) = val; end
			if n>=97
				val = str2double(S(84:97));
				if isfinite(val), tmp{6}(i,prn) = val; end
			if n>=113
				val = str2double(S(100:113));
				if isfinite(val), tmp{7}(i,prn) = val; end
			if n>=129
				val = str2double(S(116:129));
				if isfinite(val), tmp{8}(i,prn) = val; end
			if n>=145
				val = str2double(S(132:145));
				if isfinite(val), tmp{9}(i,prn) = val; end
			if n>=161
				val = str2double(S(148:161));
				if isfinite(val), tmp{10}(i,prn) = val; end
			if n>=177
				val = str2double(S(164:177));
				if isfinite(val), tmp{11}(i,prn) = val; end
			if n>=193
				val = str2double(S(180:193));
				if isfinite(val), tmp{12}(i,prn) = val; end
			if n>=209
				val = str2double(S(196:209));
				if isfinite(val), tmp{13}(i,prn) = val; end
			if n>=225
				val = str2double(S(212:225));
				if isfinite(val), tmp{14}(i,prn) = val; end
			if n>=241
				val = str2double(S(228:241));
				if isfinite(val), tmp{15}(i,prn) = val; end
			if n>=257
				val = str2double(S(244:257));
				if isfinite(val), tmp{16}(i,prn) = val; end
			end, end, end, end, end, end, end, end, end, end
			end, end, end, end, end, end, end
		end
	end
	fclose(fid);
end

% Eliminate outliers in code observables, trim out the zero tails and compact the observation arrays
for k=1:ntype
	switch GNSS
	case 'R'
		switch Type(k,:)
		case {'C1' 'C2' 'C7' 'P1' 'P2' 'C1C' 'C1P' 'C2C' 'C2P' 'C3I' 'C3Q' 'C3X'}
			tmp{k}(tmp{k}<bnd_pr.r(1) | tmp{k}>bnd_pr.r(2)) = 0;
		end
	case 'E'
		switch Type(k,:)
		case {'C1' 'C5' 'C6' 'C7' 'C8' 'C1A' 'C1B' 'C1C' 'C1X' 'C1Z' 'C5I' 'C5Q' 'C5X' 'C6A' 'C6B' 'C6C' 'C6X' 'C6Z' 'C7I' 'C7Q' 'C7X' 'C8I' 'C8Q' 'C8X'}
			tmp{k}(tmp{k}<bnd_pr.e(1) | tmp{k}>bnd_pr.e(2)) = 0;
		end
	case 'C'
		switch Type(k,:)
		case {'C1' 'C2' 'C6' 'C7' 'C1I' 'C1Q' 'C1X' 'C2I' 'C2Q' 'C2X' 'C6I' 'C6Q' 'C6X' 'C7I' 'C7Q' 'C7X'}
			tmp{k}(tmp{k}<bnd_pr.c(1) | tmp{k}>bnd_pr.c(2)) = 0;
		end
	end
	eval([Type(k,:) ' = sparse(tmp{k}(1:i,:));']);
end

% Convert into fewer-bit types
year       = int16(year'); % 2- or 4-digit year
month      = int8(month');
day        = int8(day');
hour       = int8(hour');
minute     = int8(minute');
second     = single(second');
event_flag = logical(event_flag');

% Save to file
OUTFILE = [INFILE(1:8) '.mat'];
switch ntype
case  1, save(OUTFILE,Type(1,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  2, save(OUTFILE,Type(1,:),Type(2,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  3, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  4, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  5, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  6, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  7, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  8, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case  9, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 10, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 11, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 12, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 13, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 14, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 15, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 16, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 17, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 18, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 19, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 20, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 21, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),Type(21,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 22, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),Type(21,:),Type(22,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 23, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),Type(21,:),Type(22,:),Type(23,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 24, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),Type(21,:),Type(22,:),Type(23,:),...
	Type(24,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 25, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),Type(21,:),Type(22,:),Type(23,:),...
	Type(24,:),Type(25,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
case 26, save(OUTFILE,Type(1,:),Type(2,:),Type(3,:),Type(4,:),Type(5,:),...
	Type(6,:),Type(7,:),Type(8,:),Type(9,:),Type(10,:),Type(11,:),...
	Type(12,:),Type(13,:),Type(14,:),Type(15,:),Type(16,:),Type(17,:),...
	Type(18,:),Type(19,:),Type(20,:),Type(21,:),Type(22,:),Type(23,:),...
	Type(24,:),Type(25,:),Type(26,:),...
	'year','month','day','hour','minute','second','event_flag','GNSS',...
	'rinex_version','Domes','x','y','z','intv','TimeSystem');
end
if rv==3
	if exist('receiver_clk_offset','var')
		save('-append',OUTFILE,'receiver_clk_offset');
	end
	if exist('sid','var')
		save('-append',OUTFILE,'sid','freqn');
	end
end
sflag = false;
disp(['FILE CREATION: ' OUTFILE]);

