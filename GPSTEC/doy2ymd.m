function [year,month,day] = doy2ymd(year,doy)
%DOY2YMD Convert day of year into month and day of the month.
% [Y,M,D] = DOY2YMD(YR,DOY) converts two-digit or four-digit year YR and the
% day of year DOY into four-digit year Y, month M, and day of the month D.
% The year is valid between 1951-2050.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2008/04/10
%   VERSION	:	1.0 2015/12/23

narginchk(2,2);

% Convert 2-digit year into 4-digit year (valid for 1951-2050)
l = year>=51 & year<100;
if any(l)
	year(l) = year(l)+1900;
end
l = year<=50;
if any(l)
	year(l) = year(l)+2000;
end

% Convert DOY into month and day
[year,month,day] = datevec(datenum(year,1,1) + doy-1);
end