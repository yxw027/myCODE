function recdcb = getrecdcb_r(Site,yr,doy,stec,el,ep,xi,yi,zi,sid,freqn,Domes)
%GETRECDCB_R   Estimate differential code bias of GLONASS receiver.
% DCB = GETRECDCB_R('SITE',Y,D,STEC,EL,EP,SID,FREQN,X,Y,Z) estimates the
% differential code biases (DCB) in nanoseconds for the GLONASS site, SITE in
% 4-char, with the given year Y and day-of-year D, the slant TEC, STEC in TECU,
% the corresponding elevation angles of satellites, EL in degrees, epochs EP in
% hours, XYZ positions of satellites in km.
%
% DCB = GETRECDCB_R(..., 'DOMES') further assigns the DOMES number.  The DOMES
% number should be in nnnnnMnnnXX or nnnnnSnnnXX format; otherwise, it will be
% discarded.

%   AUTHOR	:	Min-Yang Chou, Ho-Fang Tsai
%   SINCE	:	2015/12/18
%   VERSION	:	0.9 2016/12/22

error(nargchk(11,12,nargin));

% Load mapping package if not loaded
if exist('OCTAVE_VERSION','builtin')
	Package = pkg('list','mapping');
	if ~Package{1}.loaded
		pkg load mapping;
	end
end

% Parameters
clight = 2.99792458d8; % Light speed in m/s
f1o = 1.602d9; % Base frequency of L1 in Hz
f1i = 0.5625d6; % Frequency increment of L1 in Hz
f2o = 1.246d9; % Base frequency of L2 in Hz
f2i = 0.4375d6; % Frequency increment of L2 in Hz
ionfac = 40.3082d16; % Ionospheric factor in MKS units (d = 40.3TEC/(f^2))
rad = 180. /pi;
nsat = 26; % Number of satellites
km2m = 1e3; % Conversion factor for km to m
order = 4;

% Deal with time information
Y = num2str(yr,'%04u'); Y2 = Y(3:4);
D = num2str(doy,'%03u'); % Day of year in 3-char

% Reconstruct frequencies
k = zeros(1,nsat);
for i=1:nsat
	f = find(sid==i);
	if f
		k(i) = freqn(f(1));
	end
end
f1 = f1o+f1i*k; % L1 frequencies in Hz
f2 = f2o+f2i*k; % L2 frequencies in Hz
f12 = f1.^2; f22 = f2.^2;
m2tecu = f12.*f22./(f12-f22)/ionfac; % Conversion factor from meters to TECU
tecu2ns = 1e9/clight./m2tecu; % Conversion factor from TECU to nanoseconds

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
sg = size(glon);
A = azimuth(repmat(rlat,sg(1),sg(2)),repmat(rlon,sg(1),sg(2)),glat,glon);
%A = azimuth(rlat,rlon,glat,glon);
EL = full(el);
Idx = ~stec;
stec(Idx) = NaN;
A(Idx) = NaN;
EL(Idx) = NaN;
t_r = ep*2*pi/24;

[H,L] = Hmatrix(stec.*repmat(tecu2ns,size(stec,1),1),A,EL,t_r,rlon,rlat,order);
R = H\L;
recdcb = R(1);

