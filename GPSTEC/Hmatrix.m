function [H,L]=Hmatrix(stec,A,EL,t_r,rlon,rlat,order)
%HMATRIX  Generate the design matrix for spherical harmonics for GNSS.
% [H, L] = HMATRIX(STEC,A,EL,OS,RLON,RLAT,ORDER) processes the slant TEC data,
% STEC, azimuth and elevation angle of the ionospheric point (A and EL) in
% degrees, longitude shift OS from Universal time to local time in radians,
% longitude and latitude of the specified site (RLON, RLAT) in degrees by the
% set of the spherical harmonic functions with the give order, and output the
% design matrix H and the vertical TEC vector L in the same unit as STEC based
% on the least squares method.
%  
% NOTICE: For one-day worth of GPS/GLONASS/Galileo/BDS data only.

%   AUTHOR	:	Min-Yang Chou, Ho-Fang Tsai
%   SINCE	:	2015/12/15
%   VERSION	:	0.2 2016/08/30

% Parameters
dhour = 12; % Number of every 2-hour sets in a day
rearth = 6371; % Mean radius of earth in km
h = 325; % Height of ionosphere in km
rad = 180. /pi;

%
nsat = size(stec,2); % Number of satellites
datanumber = sum(~isnan(stec(:)));
[MAPF,zp] = slant2(rearth,90-EL,h);
ltr = length(t_r);% Number of epoches (86400 for 1-s or 2880 for 30-s sampling in a day)
dltr = ltr/dhour; % Number of 2-hour epoches
t_r = repmat(t_r,1,nsat);
H = zeros(datanumber,(order+1)^2*12+1);
L = zeros(datanumber,1);

sbm = repmat(rlat/rad,ltr,nsat);
slm = repmat(rlon/rad,ltr,nsat);
[b,s] = Get_IPP(EL/rad,A/rad,sbm,slm,zp,t_r);
c = 1;
for i = 1: dhour
	k = (i-1)*dltr+1:dltr*i; % Epoch number
	zp1 = zp(k,:);
	MAPF1 = MAPF(k,:);
	stec1 = stec(k,:);
	b1 = b(k,:);
	s1 = s(k,:);
	II = ~isnan(b1);
	b1 = b1(II);
	s1 = s1(II);
	stec1 = stec1(II);
	zp1 = zp1(II);
	MAPF1 = MAPF1(II);
	H(c:c+length(zp1)-1,1) = -MAPF1;
	st = (order+1)^2*(i-1)+2;
	ed = (order+1)^2*i+1;
	H(c:c+length(zp1)-1,st:ed) = Get_coef(b1,s1,order);
	L(c:c+length(zp1)-1,1) = stec1.*MAPF1;
	c = c+length(zp1);
end

