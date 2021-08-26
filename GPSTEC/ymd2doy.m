function doy = ymd2doy(year,month,day)
%YMD2DOY Convert year, month and day into day of year.
% DOY = YMD2DOY(Y,M,D) converts two- or four-digit year Y, month M and day of
% the month D into day of year DOY.  The year is valid between 1951-2050.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/25
%   VERSION	:	0.5 2015/12/23

error(nargchk(3,3,nargin));

% Convert 2-digit year into 4-digit year (valid for 1951-2050)
l = year>=51 & year<100;
if any(l)
	year(l) = year(l)+1900;
end
l = year<=50;
if any(l)
	year(l) = year(l)+2000;
end

% Compute day of the year
doy = datenum(year,month,day)-datenum(year,1,1)+1;

