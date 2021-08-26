function [w,wd] = doy2gps(year,doy)
%DOY2GPS Convert day of year into GPS week and day of the week.
% [W,WD] = DOY2GPS(Y,DOY) converts two-digit or four-digit year Y and the day
% of year DOY into GPS week W and day of week WD.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/05/07
%   VERSION	:	0.5 2013/12/23

error(nargchk(2,2,nargin));

% Parameters
d0 = datenum(1980,1,6,0,0,0); % GPS start time in days

% Convert into [year, month, day]
%[year,month,day] = doy2ymd(year,doy);

% Convert DOY to GPS week and day of week
gpsd = datenum(year,1,1)+doy-1-d0; % GPS days
w = floor(gpsd/7);
wd = gpsd-w*7;

