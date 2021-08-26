function [recdcb,satedcb] = getalldcb(Site,stec,el,ep,xi,yi,zi,Domes)
%GETALLDCB   Estimate differential code bias of GPS receiver and satellites.
% [RDCB, SDCB] = GETALLDCB('SITE',STEC,EL,EP,X,Y,Z) estimates the differential
% code biases (DCB) of the GPS receiver, RDCB, and satellites, SDCB, in nano-
% seconds for site, SITE in 4-char, with the given slant TEC, STEC in TECU, the
% corresponding elevation angles of satellites, EL in degrees, epochs EP in
% hours, XYZ positions of satellites in km.
%
% DCB = GETALLDCB(..., 'DOMES') further assigns the DOMES number.  The DOMES
% number should be in nnnnnMnnnXX or nnnnnSnnnXX format; otherwise, it will be
% discarded.

%   AUTHOR	:	Min-Yang Chou, Ho-Fang Tsai
%   SINCE	:	2016/11/29
%   VERSION	:	0.2 2016/12/22

error(nargchk(7,8,nargin));

% Load mapping package if not loaded
if exist('OCTAVE_VERSION','builtin')
	Package = pkg('list','mapping');
	if ~Package{1}.loaded
		pkg load mapping;
	end
end

% Parameters
clight = 2.99792458d8; % Light speed in m/s
rearth = 6371; % Mean radius of earth in km
f1 = 1.57542d9; % L1 frequency in Hz
f2 = 1.2276d9; % L2 frequency in Hz
f12 = f1^2; f22 = f2^2;
ionfac = 40.3082d16; % Ionospheric factor in m/TECU/s^2
m2tecu = f12*f22/(f12-f22)/ionfac; % Conversion factor from meters to TECU
tecu2ns = 1e9/clight/m2tecu; % Conversion factor from TECU to nanoseconds
rad = 180. /pi;
nsat = 32; % Number of satellites
km2m = 1e3; % Conversion factor for km to m
order = 4;
minobs = 100; % Minimum number of observables for a satellite

% Get the position of the specified site
[rx,ry,rz] = readcrd(Site,Domes);
[rlat,rlon] = ecef2lla(rx*km2m,ry*km2m,rz*km2m);
rlat = rlat*rad;
rlon = rlon*rad;
if rlon>180
	rlon = rlon-360;
end

% Estimate the azimuth of the satellites
[glat,glon] = ecef2lla(xi*km2m,yi*km2m,zi*km2m);
glat = glat*rad;
glon = glon*rad;
I = glon>180;
glon(I) = glon(I)-360;
A = azimuth(rlat,rlon,glat,glon);
EL = full(el);
Idx = ~stec;
stec(Idx) = NaN;
A(Idx) = NaN;
EL(Idx) = NaN;
t_r = ep*2*pi/24;
II = isnan(stec);
stec1 = stec;
stec1(II) = 0;
nobs = sum(logical(stec1));
l = nobs<minobs;
stec1(:,l) = 0;
stec(:,l) = NaN;
A(:,l) = NaN;
EL(:,l) = NaN;
satidx = find(sum(stec1)); % Observing satellites
nsatidx = length(satidx); % Number of observing satellites
[H,L] = Hmatrix_m(stec*tecu2ns,A,EL,t_r,rlon,rlat,order);
H1 = linspace(0,0,(order+1)^2*12+nsatidx+1);
H1(1,2:nsatidx+1) = 1;
H = [H; H1];
L =[L; 0];
R = H\L;
recdcb = R(1);
satedcb1 = R(2:nsatidx+1);
satedcb = zeros(1,nsat); % Use zeros() instead of nan()
satedcb(1,satidx) = satedcb1;
