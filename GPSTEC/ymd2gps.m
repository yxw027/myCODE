function [w,wd] = ymd2gps(year,month,day)
%YMD2GPS   Convert year, month and day into GPS week and day of the week.
% [W,WD] = YMD2GPS(Y,M,D) converts two- or four-digit year Y, month M and day
% of the month D into GPS week W and day of week WD.  The year is valid between
% 1951-2050.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/12
%   VERSION	:	0.5 2015/12/23

narginchk(3,3);

% Parameters
d0 = datenum(1980,1,6,0,0,0); % GPS start time in days

% Convert 2-digit year into 4-digit year (valid for 1951-2050)
l = year>=51 & year<100;
if any(l)
	year(l) = year(l)+1900;
end
l = year<=50;
if any(l)
	year(l) = year(l)+2000;
end

% Conversion
gpsd = datenum(year,month,day)-d0; % GPS days
w = floor(gpsd/7);
wd = gpsd-w*7;
end
