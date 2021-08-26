function [H,L]=Hmatrix_m(stec,A,EL,t_r,rlon,rlat,order)
%HMATRIX_M  Generate the design matrix for spherical harmonics for GPS and BDS.
% [H, L] = HMATRIX_M(STEC,A,EL,OS,RLON,RLAT,ORDER) processes the slant TEC
% data, STEC, azimuth and elevation angle of the ionospheric point (A and EL)
% in degrees, longitude shift OS from Universal time to local time in radians,
% longitude and latitude of the specified site (RLON, RLAT) in degrees by the
% set of the spherical harmonic functions with the give order, and output the
% design matrix H and the vertical TEC vector L in the same unit as STEC based
% on the least squares method.
%  
% NOTICE: For one-day worth of BDS data only.

%   AUTHOR	:	Min-Yang Chou, Ho-Fang Tsai
%   SINCE	:	2016/06/24
%   VERSION	:	0.4 2016/12/22

% Parameters
dhour = 12; % Number of every 2-hour sets in a day
rearth = 6371; % Mean radius of earth in km
h = 450; % Height of ionosphere in km
rad = 180. /pi;

%
nsat = size(stec,2); % Number of satellites
datanumber = sum(~isnan(stec(:)));
[MAPF,zp] = slant2(rearth,90-EL,h);

ltr = length(t_r);% Number of epoches (86400 for 1-s or 2880 for 30-s sampling in a day)
dltr = ltr/dhour; % Number of 2-hour epoches
t_r = repmat(t_r,1,nsat);
II = isnan(stec);
stec1 = stec;
stec1(II) = 0;
satidx = find(sum(stec1)); % Observing satellite id
nosatidx = find(~sum(stec1)); % Non-observing satellite id
nsatidx = size(satidx,2);
H = zeros(datanumber,(order+1)^2*12+1+nsatidx);
L = zeros(datanumber,1);
satidx11 = repmat(1:nsatidx,ltr,1);% Satellite index position
satidx1 = zeros(size(stec,1),size(stec,2));
for i =1:nsatidx
    satidx1(:,satidx(i)) = i;
end
sbm = repmat(rlat/rad,ltr,nsat);
slm = repmat(rlon/rad,ltr,nsat);
[b,s] = Get_IPP(EL/rad,A/rad,sbm,slm,zp,t_r);
c = 1;
kk = 1;
for i=1:dhour
	k = (i-1)*dltr+1:dltr*i; % Epoch number
	MAPF1 = MAPF(k,:);
	stec1 = stec(k,:);
	b1 = b(k,:);
	s1 = s(k,:);
        satidx2 = satidx1(k,:);
	II = ~isnan(b1);
	b1 = b1(II);
	s1 = s1(II);
	stec1 = stec1(II);
	MAPF1 = MAPF1(II);
        satidx2 = satidx2(II);
	H(c:c+length(stec1)-1,1) = -MAPF1;
        usat2 = unique(satidx2);
        cc =1;
        for k2 = 1:length(usat2)
            lidx = sum(satidx2 ==usat2(k2)); 
            H(kk:kk+lidx-1,usat2(k2)+1) = -MAPF1(cc:cc+lidx-1);
            kk = kk+lidx;
            cc = cc+lidx;
        end 
	st = (order+1)^2*(i-1)+2+nsatidx;
	ed = (order+1)^2*i+1+nsatidx;
	H(c:c+length(stec1)-1,st:ed) = Get_coef(b1,s1,order);
	L(c:c+length(stec1)-1,1) = stec1.*MAPF1;
	c = c+length(stec1);
end

