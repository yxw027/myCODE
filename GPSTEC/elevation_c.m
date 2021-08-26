function el = elevation_c(Site,x,y,z,minel,Domes)
%ELEVATION_C   Calculate elevation angle for BDS.
% EL = ELEVATION_C('SITE',X,Y,Z,MINEL) calculates the elevation angles EL of
% the satellites at [X, Y, Z] with respect to BDS site SITE in 4-char above the
% minimum (cutoff) elevation angle MINEL.  All angles are in degrees and all
% lengths are in km, ECF coordinate.
%
% For example, ELEVATION_C('coco',X,Y,Z,20) returns elevation angles greater
% than 20 degrees for site 'coco'.  Any elevation angles less than MINEL are
% set to be zero.
%
% EL = ELEVATION_E('SITE',X,Y,Z,MINEL,'DOMES') further assigns the DOMES
% number.  The DOMES number should be in nnnnnMnnnXX or nnnnnSnnnXX format;
% otherwise, it will be discarded.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2016/06/24
%   VERSION	:	0.4 2016/12/09

error(nargchk(5,6,nargin));

% Parameters
rad = 180. /pi;
nsat = 35; % Number of BDS satellites

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

% Transport arrays if not (n x 35)
[nr,nc] = size(x);
if nr==nsat
	x = x'; y = y'; z = z';

	% Swap the numbers of rows and columns
	n = nr;
	nr = nc;
	nc = n;
elseif nc~=nsat
	error('XYZ should be in (n x 35) or (35 x n).');
end

% Get the position vector of the specified site
[rx,ry,rz] = readcrd(Site,Domes); % XYZ in km
r = [rx ry rz]; % Receiver vector in km
rn = norm(r); % Length of vector 'r'
r1 = ones(nr,1)*r/rn; % Unit vector of 'r'
r = ones(nr,1)*r;

% Allocate space for new variables
el = zeros(nr,nc);

% Estimate the elevation angle
for i=1:nc
	s = [x(:,i) y(:,i) z(:,i)]; % Satellite vectors in km
	sr = s-r; % Vector of satellite-receiver links in km
	sr1 = sr./(sqrt(sr(:,1).^2+sr(:,2).^2+sr(:,3).^2)*ones(1,3)); % Unit vector of 'sr'
	el(:,i) = asin(dot(r1,sr1,2)); % Elevation angle in radians
end

% Post-process 'el'
el = el*rad; % Convert radians into degrees
el(isnan(el)) = 0; % Set all NaN's to be zeros
el(el<minel) = 0; % Set all elevation angles less then MINEL to be zeros
%el = sparse(el); % Remove zeros
