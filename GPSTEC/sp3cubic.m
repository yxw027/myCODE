function sp3cubic(INFILE,rate)
%SP3CUBIC  Apply cubic interpolation on SP3 data.
% SP3CUBIC('INFILE',R) interpolates the SP3 data from 15-minute time interval
% to R second(s) sampling with cubic interpolation.  For example,
% SP3CUBIC('igs16265.mat',1) retrieves XYZ coordinates from 'igs16265.mat',
% interpolates into 1-second sampling and then writes to 'igs16265i01.mat'.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/12
%   VERSION	:	1.5 2016/06/16

error(nargchk(2,2,nargin));

% Parameters
R = num2str(rate,'%02u');
 
% Load SP3 file in MAT format
load(INFILE,'years','months','days','hours','minutes','seconds','x','y','z');

% Cubic interpolation
t = datenum(double(years),double(months),double(days),double(hours),...
	double(minutes),double(seconds))*86400; % Original epoches in seconds
[theta,phi,r] = cart2sph(x,y,z); % Theta is related to longitude (but not equal); phi is related to latitude (but not equal)
ti = (t(1):rate:t(end)+t(2)-t(1)-rate)'; % Interpolated epoches in seconds
sti = size(ti);
for i=1:size(x,2)

	% Connect any angle jump in theta
	f = find(diff(theta(:,i))<-pi); % Find any jumps less than -pi
	for j=f'
		theta(j+1:end,i) = theta(j+1:end,i)+2*pi;
	end
	f = find(diff(theta(:,i))>pi); % Find any jumps greater than pi
	for j=f'
		theta(j+1:end,i) = theta(j+1:end,i)-2*pi;
	end

	% Cubic interpolation for (theta, phi, r)
	if sum(isfinite(theta(:,i)))>2
		thetai(:,i) = interp1(t,theta(:,i),ti,'pchip','extrap'); % Interpolated theta in radians
		phii(:,i)   = interp1(t,phi(:,i),ti,'pchip','extrap'); % Interpolated phi in radians
		ri(:,i)     = interp1(t,r(:,i),ti,'pchip','extrap'); % Interpolated r
	else
		thetai(:,i) = NaN(sti);
		phii(:,i)   = thetai(:,i);
		ri(:,i)     = thetai(:,i);
	end
end
[xi,yi,zi] = sph2cart(thetai,phii,ri); % Interpolated XYZ
[years,months,days,hours,minutes,seconds] = datevec(ti/86400);
 
% Convert doubler-precision arrays into integer arrays
years   = int16(years);
months  = int8(months);
days    = int8(days);
hours   = int8(hours);
minutes = int8(minutes);
seconds = int8(seconds);

% Write to file
OUTFILE = [INFILE(1:end-4) 'i' R '.mat'];
save(OUTFILE,'xi','yi','zi','years','months','days','hours','minutes',...
	'seconds');
disp(['FILE CREATION: ' OUTFILE]);

