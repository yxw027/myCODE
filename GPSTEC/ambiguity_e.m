function [ep,stecl,amb] = ambiguity_e(Site,yr,doy,hr,el,recdcb,galdcb,dflag,rate)
%AMBIGUITY_E   Solve ambiguity for Galileo.
% [EP,STEC,AMB] = AMBIGUITY_E('SITE',Y,D,H,EL,RDCB,SDCB,F,R) solves the
% ambiguity in R-second sampling for Galileo with receiver's and satellite's
% differential code biases (RDCB and SDCB) in nanoseconds, elevation angle EL
% in degrees, for site SITE in 4-char, year Y, day-of-year D and hour H, and
% outputs the ambiguity AMB in TECU and the ambiguity-free slant TEC (STEC) in
% TECU for epoches EP in fractional hours.
%
% [EP,STEC,AMB] = AMBIGUITY_E('SITE',Y,D,[H1 H2],EL,RDCB,SDCB,F,R) specifies
% the period in hours [H1 H2].
%
% F = 1 indicates daily solution; otherwise, hourly solution.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/11/29
%   VERSION	:	2.9 2016/12/10

error(nargchk(9,9,nargin));

% Parameters
f1 = 1.57542d9; % L1 frequency in Hz
f5 = 1.17645d9; % L5 frequency in Hz
f12 = f1^2; f52 = f5^2;
ionfac = 40.3082d16; % Ionospheric factor in m/TECU/s^2
m2tecu = f12*f52/(f12-f52)/ionfac; % Conversion factor from meters to TECU
clight = 2.99792458d8; % Light speed in m/s
ns2m = clight*1e-9; % Conversion factor from nanoseconds to meters
lambda1 = clight/f1; % L1 wavelength in m
lambda5 = clight/f5; % L5 wavelength in m
re = 6371; % Mean radius of the earth in km
h1 = 650; % Maximum height of ionosphere in km
h2 = 250; % Minimum height of ionosphere in km
rad = 180. /pi;
nepoch = 3600/rate; % Number of epochs in an hour
nsat = 30; % Number of satellites
%minstec = -100; % Minimum tolerance of slant TEC in TECU
%maxstec = 500; % Maximum tolerance of slant TEC in TECU

% Set the maximum tolerance of TEC difference in TEC unit between two adjacent epoches
if rate==1
	maxtecdiff = 0.2;
else
	maxtecdiff = 1.2;
end

% Deal with time information
Y = num2str(yr,'%04u'); Y2 = Y(3:4);
D = num2str(doy,'%03u');
[year,month,day] = doy2ymd(yr,doy);
M = num2str(month,'%02u');
if size(hr,2)==1
	hr = [hr hr+1];
end
d0 = datenum(year,month,day,0,0,0); % Initial day number in days

% Initialize variables
l1 = []; l5 = l1; c1 = l1; c5 = l1; epo = l1;

% Load data
if dflag==1 % Daily
	INFILES = [Site D '0.mat']; % SSSSDDD0.mat
elseif dflag==0 % Hourly
	for i=hr(1):hr(2)-1
		TMPFILE = [Site D char(97+i) '.mat']; % SSSSDDDA.mat
		if exist(TMPFILE,'file')
			INFILES(i-hr(1)+1,:) = TMPFILE; % Collect valid file names
		end
	end
end
for i=1:size(INFILES,1)
	load(INFILES(i,:),'rinex_version');
	rv = floor(rinex_version);
	switch rv
	case 2
		load(INFILES(i,:),'C1','L1','C5','L5','year','month','day',...
			'hour','minute','second','event_flag');
		l1 = [l1;L1(~event_flag,:)];
		l5 = [l5;L5(~event_flag,:)];
		c1 = [c1;C1(~event_flag,:)];
		c5 = [c5;C5(~event_flag,:)];
		epo = [epo;datenum(double(year(~event_flag))+2e3,double(month(~event_flag)),double(day(~event_flag)),double(hour(~event_flag)),double(minute(~event_flag)),double(second(~event_flag)))]; % Epoches in days
		clear C[15] L[15] year month day hour minute second event_flag;
	case 3
		load(INFILES(i,:),'C1C','L1C','C5Q','L5Q','year','month',...
			'day','hour','minute','second','event_flag');
		l1 = [l1;L1C(~event_flag,:)];
		l5 = [l5;L5Q(~event_flag,:)];
		c1 = [c1;C1C(~event_flag,:)];
		c5 = [c5;C5Q(~event_flag,:)];
		epo = [epo;datenum(double(year(~event_flag)),double(month(~event_flag)),double(day(~event_flag)),double(hour(~event_flag)),double(minute(~event_flag)),double(second(~event_flag)))]; % Epoches in days
		clear C1C C5Q L1C L5Q year month day hour minute second event_flag;
	otherwise
		error('The RINEX version of the input file should be 2.00 or later.');
	end
end

% Allocate space for observational arrays
ep = (hr(1):rate/3600:hr(2)-rate/3600)'; % New epochs in hours
nep = size(ep,1);
L1 = zeros(nep,nsat); L5 = L1; C1 = L1; C5 = L1;

% Pre-process the observation data
epo = (epo-d0)*24; % Original epoches in hours
n = round((epo-hr(1))*nepoch)+1; % Indices for 1-s or 30-s sampling
L1(n,:) = full(l1); L1(~L1) = NaN;
L5(n,:) = full(l5); L5(~L5) = NaN;
C1(n,:) = full(c1); C1(~C1) = NaN;
C5(n,:) = full(c5); C5(~C5) = NaN;
clear c[15] l[15];

% Remove the data of invalid elevation angles
if any(any(el))
	el = full(el);
	l = ~el;
	L1(l) = NaN; L5(l) = NaN; C1(l) = NaN; C5(l) = NaN; el(l) = NaN;
else
	error('All zeros in variable el.');
end

% Deal with satellite DCBs
nr = size(C5,1);
galdcb = ones(nr,1)*galdcb*ns2m; % Convert nanoseconds into meters

% Estimate code and phase TEC respectively
recdcb = recdcb*ns2m; % Convert nanoseconds into meters
stecp = (C5-C1+galdcb+recdcb)*m2tecu; % Slant TEC based on code observables
stecl = (lambda1*L1-lambda5*L5)*m2tecu;

% Remove outliers from (code) slant TEC
%stecp(stecp<minstec | stecp>maxstec) = NaN;

% Remove isolated data (no data before and after the point(s))
for i=1:nsat
	for j=2:nr-1
		if isfinite(stecl(j-1:j+1,i))==[0 1 0]'
			stecl(j,i) = NaN;
		end
	end
	for j=2:nr-2
		if isfinite(stecl(j-1:j+2,i))==[0 1 1 0]'
			stecl(j:j+1,i) = NaN;
		end
	end
end

% Fill temporary data into small data gaps to link the broken data (will remove it after solving the ambiguity)
cflag = false(nr,nsat); % Flag array for data gaps to be filled
for i=1:nsat
	for j=2:nr-2
		if isnan(stecl(j-1:j+2,i))==[0 1 0 0]'
			stecl(j,i) = stecl(j+1,i)*2-stecl(j+2,i);
			stecp(j,i) = (stecp(j-1,i)+stecp(j+1,i))/2;
			cflag(j,i) = true;
		end
	end
	for j=3:nr-1
		if isnan(stecl(j-2:j+1,i))==[0 0 1 0]'
			stecl(j,i) = stecl(j-1,i)*2-stecl(j-2,i);
			stecp(j,i) = (stecp(j-1,i)+stecp(j+1,i))/2;
			cflag(j,i) = true;
		end
	end
	for j=3:nr-2
		if isnan(stecl(j-2:j+2,i))==[0 0 1 1 0]'
			stecl(j,i) = stecl(j-1,i)*2-stecl(j-2,i);
			stecl(j+1,i) = stecl(j+2,i);
			stecp(j,i) = (stecp(j-1,i)+stecp(j+2,i))/2;
			stecp(j+1,i) = stecp(j,i);
			cflag(j:j+1,i) = true;
		end
	end
end

% Detect and offset cycle slips
dstecl = diff(stecl);
S = slant(re,el,h1,h2);
for i=1:nsat
	for j=1:nr-1
		if abs(dstecl(j,i))>maxtecdiff*S(j,i) % Treat large TEC jump as a cycle slip
			stecl(j+1:end,i) = stecl(j+1:end,i)-dstecl(j,i); % Offset slant TEC behind this point with the TEC difference between this point and the next point
		end
	end
end

amb = stecp-stecl; % Ambiguity in TECU

for i=1:nsat
	fa = find(isfinite(amb(:,i)));
	dfa = diff(fa);
	fdfa = find(dfa>1);
	for j=0:size(fdfa,1)
		if ~size(fdfa,1)
			k = fa;
		elseif j==0
			k = fa(1:fdfa(1));
		elseif j==size(fdfa,1)
			k = fa(fdfa(end)+1:end);
		else
			k = fa(fdfa(j)+1:fdfa(j+1));
		end
		amb(k,i) = mean(amb(k,i));
	end
end

% Convert NaN's into zeros
l = isnan(stecl) | isnan(amb) | cflag;
stecl(l) = 0;
amb(l) = 0;

stec = stecl+amb;

