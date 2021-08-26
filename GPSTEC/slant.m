function s=slant(r,el,h1,h2)
%SLANT  Slant function.
% S = SLANT(R,EL,H1,H2) converts mean radius R in km, elevation angle EL in
% degrees and max/min ionospheric heights H1 and H2 in km into slant factor S,
% where vertical observable = vertical observable / S.  For ground-based
% receivers, R = 6371 km, H1 = 450~650 km, H2 = 200~250 km; for LEO recievers,
% R = 6371+h, H1 = maximum plasmaspheric height-h, H2 = h-h =0,  where h is
% altitude of the LEO orbit; H1 and H2 are counted from this altitude.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2009/02/09
%   VERSION	:	0.9 2015/12/02

error(nargchk(4,4,nargin));

% parameters
rad = 180. /pi;

r2 = r.^2; % r square
r2sin2e = r2.*sin(el/rad).^2;
s = (sqrt(r2sin2e-r2+(r+h1).^2)-sqrt(r2sin2e-r2+(r+h2).^2))./(h1-h2);
