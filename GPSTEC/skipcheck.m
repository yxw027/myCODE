function sflag = skipcheck(SiteName,yr,doy)
%SKIPCHECK   Check if the station is in the skip list.
% FLAG = SKIPCHECK('SITENAME',Y,D) checks if the input station and date are in
% the skip list STATION.CRX.  FLAG is true if so; otherwise, false.
%
% SITENAME must be in 'SSSS nnnnnMnnnXX' or 'SSSS nnnnnSnnnXX' format, where
% SSSS is the station name (case-insensitive), following by the DOMES number.
% 'XX' and even the whole DOMES number can be neglected if not available.
% Y is two-digit or four-digit year; D is the day-of-year.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2015/01/23
%   VERSION	:	0.3 2015/08/28

error(nargchk(3,3,nargin));

% Parameters
sflag = false;

% Compensate the input site name with whitespaces if necessary
SiteName = upper(SiteName);
ln = size(SiteName,2);
if ln<16
	SiteName(ln+1:16) = char(32);
end
[year,month,day] = doy2ymd(yr,doy);
date_n = datenum([year month day]);

% Load the station problem file
INFILE = 'STATION.CRX';
[fid, Msg] = fopen(INFILE,'r');
if fid==-1
	error(Msg);
end

% Read in the header of the file
for i=1:5
	fgetl(fid);
end

% Read in the body of the file
while ~feof(fid)
	S = fgetl(fid);
	if size(S,2)<54
		break;
	end
	ReceiverName = S(1:16);
	drng = [datenum([S(29:30) '/' S(32:33) '/' S(24:27)]) ...
		datenum([S(50:51) '/' S(53:54) '/' S(45:48)])]; % Date range
	if strcmp(SiteName,ReceiverName) && date_n>=drng(1) && date_n<=drng(2)
		sflag = true;
		break;
	end
end
fclose(fid);

