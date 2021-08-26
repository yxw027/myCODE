function [year,month,day] = gps2ymd(w,wd)
%GPS2YMD   Convert GPS week into year, month and day.
% [Y,M,D] = GPS2YMD(W,WD) converts the GPS week W and the day of week WD into
% year Y, month M and day of the month D.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/12
%   VERSION	:	0.5 2015/12/23

error(nargchk(2,2,nargin));

% Parameters
d0 = datenum(1980,1,6,0,0,0); % GPS start time in days

% Conversion
[year,month,day] = datevec(d0+w*7+wd);

