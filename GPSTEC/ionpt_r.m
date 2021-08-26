function [lat,lon] = ionpt_r(Site,x,y,z,el,ih,Domes)
%IONPT_R  Estimate location of the ionospheric points for GLONASS.
% [LAT,LON] = IONPT_R('SITE',X,Y,Z,EL,IH) estimates the latitutde LAT and
% longitude LON of the ionospheric points in ITRF2000 ellipsoid coordinate with
% respect to the satellites at [X, Y, Z] in km, ECF coordinate, and the GLONASS
% site,  SITE in 4-char, for the elevation angle EL in degrees and the mean
% ionospheric height IH in km.
%
% [LAT,LON] = IONPT_R('SITE',X,Y,Z,EL,IH,'DOMES') further assigns the DOMES
% number.  The DOMES number should be in nnnnnMnnnXX or nnnnnSnnnXX format;
% otherwise, it will be discarded.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2013/03/21
%   VERSION	:	1.2 2016/12/09

error(nargchk(6,7,nargin));

% Parameters
km2m = 1d3; % Factor between km and meter
rad = 180. /pi;
x = x*km2m; y = y*km2m; z = z*km2m; % Convert km into meters
ih = ih*km2m; % Convert km into meters
nsat = 26; % Number of GLONASS satellites

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

% Transport arrays if not (n x 26)
[nr,nc] = size(x);
if nr==nsat
	x = x'; y = y'; z = z';

	% Swap the numbers of rows and columns
	n = nr;
	nr = nc;
	nc = n;
elseif nc~=nsat
	error('XYZ should be in (n x 26) or (26 x n).');
end

% Get the position vector of the specified site
[rx,ry,rz] = readcrd(Site,Domes); % XYZ in km
r = [rx ry rz]*km2m; % Receiver vector in meters
rn = norm(r); % Length of vector 'r'

% Transport the elevation array if not (n x 26)
[nr,nc] = size(el);
if nr==nsat
	el = el';

	% Swap the numbers of rows and columns
	n = nr;
	nr = nc;
	nc = n;
elseif nc~=nsat
	error('EL should be in (n x 26) or (26 x n).');
end

% Calculate distance between the ionospheric point and the receiver
el = full(el)/rad; % Convert degrees into radians
minel = min(min(el(el>0))); % Minimum elevation angle in radians
arc = acos(rn*cos(el)/(rn+ih))-el; % Central angle between satellites and the receiver
dip = sqrt(rn^2+(rn+ih)^2-2. *rn*(rn+ih)*cos(arc)); % Distance between the ionospheric point and the receiver in meters

% Allocate space for new variables
xi = zeros(nr,nc);
yi = xi; zi = xi;

% Estimate the distance between satellites and the receiver
r = ones(size(x,1),1)*r;
for i=1:nc
	s = [x(:,i) y(:,i) z(:,i)]; % Satellite array in meters
	sr = s-r; % Vector of satellite-receiver links
	srn(:,i) = sqrt(sr(:,1).^2+sr(:,2).^2+sr(:,3).^2); % Length of 'sr'
end

% Calculate XYZ coordinate of the ionospheric points
ir = dip./srn;
xi = r(1,1)+ir.*(x-r(1,1)); % in meters
yi = r(1,2)+ir.*(y-r(1,2));
zi = r(1,3)+ir.*(z-r(1,3));

% Convert ECEF XYZ into ITRF2000 coordinate
[lat,lon,alt] = ecef2lla(xi,yi,zi);

% Set elements to be zero for all elevation < minimum(elev)
l = el<minel;
lat(l) = 0;
lon(l) = 0;
lat = lat*rad;
lon = lon*rad;
l = lon>180;
lon(l) = lon(l)-360; % Convert (180 E, 360 E] into [-180 E, 0)
%lon = sparse(lon);
