function [s,zp]=slant2(r,z,H);
%SLANT2  Slant function 2.
% [S, ZP] = SLANT2(R,Z,H) converts mean radius R in km, ground zenith angle Z
% in degrees and the height of the ionospheric point H in km into satellite's
% zenith angle ZP in radians and the slant factor S.  For ground-based
% receivers, the mean Earth radius R = 6371 km, and H = 450 km.

%   AUTHOR	:	Min-Yang Chou
%   SINCE	:	2015/12/15
%   VERSION	:	0.3 2016/11/08

error(nargchk(3,3,nargin));

zp = asin(r./(r+H).*sind(0.9782.*z));
s = cos(zp);
