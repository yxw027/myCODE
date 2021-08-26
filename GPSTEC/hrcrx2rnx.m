function hrcrx2rnx(Site,yr,doy,hr,rate)
%HRCRX2RNX  Convert high-rate zipped compact RINEX OBS files into RINEX files.
% HRCRX2RNX('SITE',Y,D,H,R) converts four 15-minute worth of zipped compact
% RINEX OBS files in one-second sampling into one-hour worth of RINEX files in
% R-second sampling for site SITE in 4-char, year Y, day-of-year D, and
% specified hour H.
%
% HRCRX2RNX('SITE',Y,D,[H1 H2],R) specifies the period in hours [H1 H2].
%
% The source file should be named in SSSSDDDAMM.YYd.Z, where DDD is 3-digit day
% of year, A is the hour alphabet (either uppercase or lowercase), MM is 2-
% digit minutes (00, 15, 30, or 45 only), YY is 2-digit year.  The output name
% is in form of SSSSDDDA.YYo.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2011/03/11
%   VERSION	:	2.5 2016/08/28

error(nargchk(5,5,nargin));

% Parameters
Y = num2str(yr,'%04u'); Y2 = Y(3:4);
D = num2str(doy,'%03u'); % Day of year in 3-char
R = int2str(rate); % Rate in char
if isunix || ismac
	CMD = 'grep -v ';
elseif ispc
	CMD = 'findstr /V ';
end

% Start hourly looping
for hour=hr(1):hr(end)
	fflag = false; % File flag
	H = num2str(hour,'%02u'); % 2-char hour
	A = char(hour+97); % 1-char hour alphabet (a=0; b=1; c=2, etc.)
	INPATH = pwd;
%	INPATH = ['/pub/gnss/igs/data/highrate/' Y '/' D '/' Y2 'd/' H];

	% Start 15-minute looping
	for minute=0:15:45
		M = num2str(minute,'%02u'); % Minute in 2-char
		INFILE = [Site D A M '.' Y2 'd.Z']; % Format SSSSDDDA[00|15|30|45].YYd.Z

		% Convert any zipped compact RINEX files into RINEX files
		if exist([INPATH '/' INFILE],'file')
			TMPFILE = [Site D A M '.' Y2 'o']; % Format SSSSDDDA[00|15|30|45].YYo
			system(['zcat ' INPATH '/' INFILE '|crx2rnx ->' TMPFILE]);
			fflag = true; % File exists
		end
	end

	% Reset interval and remove any comment rows
	if fflag
		OUTFILE = [Site D A '.' Y2 'o']; % SSSSDDDA.YYo

		% Combine SSSSDDDA[00|15|30|45].YYo to SSSSDDDA.YYo.tmp
		system(['teqc -phc -O.dec ' R 's ' Site D A '??.' Y2 'o > ' OUTFILE '.tmp']);

		% Remove unwanted rows
		system([CMD '"^                            4  1" ' OUTFILE '.tmp | ' CMD '"^ .*\.0000000  .  .$" | ' CMD '"COMMENT$" > ' OUTFILE]); % SSSSDDDA.YYo.tmp to SSSSDDDA.YYo
		pause(5); % Make sure grep processing complete 

		% Remove temporary files
		delete([Site D A '??.' Y2 'o']);
		delete([OUTFILE '.tmp']);

		disp(['FILE CREATION: ' OUTFILE]);
	else
		warning([INFILE 'Not found in ' INPATH '.']);
	end
end

