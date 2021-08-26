function [b,s] = Get_IPP(E,A,B,L,z,t_r)
%GET_IPP    Estimate the location of the ionospheric points.
% [LAT,LON] = GET_IPP(EL,A,B,L,Z,OS) estimates the latitutde LAT and the sun
% -fixed longitude LON of the ionospheric points with the given elevation angle
% EL in radians, azimuth A in radians, the receiver's latitude and latitude (B
% and L) in radians, the zenith angle at the ionospheric point Z, and the
% longitude shift OS in radians from universal time to local time.

%   AUTHOR	:	Min-Yang Chou
%   SINCE	:	2015/12/15
%   VERSION	:	0.1 2015/12/18

t = pi/2-E-z;
b = asin(sin(B).*cos(t)+cos(B).*sin(t).*cos(A));
s = L+asin(sin(t).*sin(A)./cos(t));
s = s+t_r-pi;

