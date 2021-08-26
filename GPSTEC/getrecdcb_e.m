function recdcb = getrecdcb_e(Site,stec,el,ep,xi,yi,zi,Domes)
%GETRECDCB_E   Estimate differential code bias of Galileo receiver.
% DCB = GETRECDCB_E('SITE',STEC,EL,EP,X,Y,Z) estimates the differential code
% biases (DCB) in nanoseconds for the Galileo site, SITE in 4-char, with the
% given slant TEC, STEC in TECU, the corresponding elevation angles of
% satellites, EL in degrees, epochs EP in hours, XYZ positions of satellites in
% km.
%
% DCB = GETRECDCB_E(..., 'DOMES') further assigns the DOMES number.  The DOMES
% number should be in nnnnnMnnnXX or nnnnnSnnnXX format; otherwise, it will be
% discarded.

%   AUTHOR	:	Min-Yang Chou, Ho-Fang Tsai
%   SINCE	:	2015/12/18
%   VERSION	:	0.8 2016/12/22

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
f1 = 1.57542d9; % L1 frequency in Hz
f5 = 1.17645d9; % L5 frequency in Hz
f12 = f1^2; f52 = f5^2;
ionfac = 40.3082d16; % Ionospheric factor in m/TECU/s^2
m2tecu = f12*f52/(f12-f52)/ionfac; % Conversion factor from meters to TECU 
tecu2ns = 1e9/clight/m2tecu; % Conversion factor from TECU to nanoseconds
rad = 180. /pi;
nsat = 30; % Number of satellites
km2m = 1e3; % Conversion factor for km to m
order = 4;

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

[H,L] = Hmatrix(stec*tecu2ns,A,EL,t_r,rlon,rlat,order);
R = H\L;
recdcb = R(1);

