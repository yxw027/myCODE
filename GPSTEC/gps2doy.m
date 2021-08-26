function [year,doy] = gps2doy(w,wd)
%GPS2xDOY   Convert GPS week into day of year.
% [Y,DOY] = GPS2YMD(W,WD) converts the GPS week W and the day of week WD into
% year Y and day of the year DOY.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2016/12/19
%   VERSION	:	0.1 2016/12/19

error(nargchk(2,2,nargin));

% Parameters
d0 = datenum(1980,1,6,0,0,0); % GPS start time in days

% Conversion
[year,month,day] = datevec(d0+w*7+wd);
doy = ymd2doy(year,month,day);

