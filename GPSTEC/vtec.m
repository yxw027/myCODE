function [v,mintec,rlat,rlon] = vtec(Site,doy,ep,stec,el,rate,ippalt,Domes)
%VTEC  Calculate vertical TEC.
% [VTEC,STEC,MINTEC,LAT,LON] = VTEC('SITE',D,EP,STEC,EL,R,ALT) converts the slant
% TEC into vertical TEC for the specified site.  The input includes:
%
% SITE: Site name in 4-char
% D: Integer day of year
% EP: Epoches in fractional hours
% STEC: Slant TEC in TEC unit
% EL: Elevation angle of satellites in degrees
% R: Sampling rate in seconds
% ALT: Ionospheric points altitude
% It returns:
%
% VTEC:  Vertical TEC in TEC unit
% MINTEC:  Minimum TEC in TEC unit
% LAT, LON:  Latitude and longitude of the receiver in degrees N and E
% 
% [...] = VTEC(...,'DOMES') further assigns the DOMES number.  The DOMES number
% should be in nnnnnMnnnXX or nnnnnSnnnXX format; otherwise, it will be
% discarded.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/12
%   VERSION	:	3.3 2016/08/09

error(nargchk(7,8,nargin));

% Parameters
D = num2str(doy,'%03u'); % Day of year in 3-char
re = 6371; % Mean radius of the earth in km
h1 = 650; % Upper bound of the ionospheric height in km
h2 = 250; % Lower bound of the ionospheric height in km
Site_lc = lower(Site); % Lower cases
rad = 180. /pi;
km2m = 1d3; % Factor between km and meter
ptmin = [4 4.5]; % Period of TEC minimum in hours

% Deal with the DOMES number
if exist('Domes','var')
	nDomes = size(Domes,2);
	if nDomes>=9 && nDomes<11
		Domes = [Domes char(32*ones(1,11-nDomes))]; % Compensate the tail with whitespaces
	else
		Domes = Domes(1,1:11); % Trim out the tail if any
	end
else
	Domes = char(32*ones(1,11)); % Assign all whitespaces to Domes
end

% Calcualate VTEC
el = full(el);
el(~el) = NaN;
%S = slant(re,el,h1,h2);
S = slant2(re,90-el,ippalt);
stec(stec<0) = 0; % Set 0 for all slant TEC < 0
v = stec.*S; % Vertical TEC in TECU

% Set VTEC values more than 200 TECU to be zeros
l = v>200 | isnan(v);
v(l) = 0;

% Get the position of the specified site
[rx,ry,rz] = readcrd(Site,Domes); % XYZ in km
[rlat,rlon] = ecef2lla(rx*km2m,ry*km2m,rz*km2m);
rlon = rlon*rad; % Convert degrees into radians
rlat = rlat*rad; % Convert degrees into radians

% Shift longitude range from (180+, 360] to (-180-,0]
if rlon>180
	rlon = rlon-360;
end

% Estimate the minimum TEC
ep_lt = ep+rlon/15; % Convert GPS time (almost UT) into local time in hours
l = (ep_lt>ptmin(1) & ep_lt<ptmin(2)) | (ep_lt>ptmin(1)-24 & ...
	ep_lt<ptmin(2)-24) | (ep_lt>ptmin(1)+24 & ep_lt<ptmin(2)+24);
v4 = v(l,:);
mintec = mean(v4(logical(v4)));
disp(['Minimum TEC = ' num2str(mintec,'%6.2f') ' TECU.']);

