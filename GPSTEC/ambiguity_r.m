function [ep,stecl,amb] = ambiguity_r(Site,yr,doy,hr,el,recdcb,glodcb,dflag,rate)
%AMBIGUITY_R   Solve ambiguity for GLONASS.
% [EP,STEC,AMB] = AMBIGUITY_R('SITE',Y,D,H,EL,RDCB,SDCB,F,R) solves the
% ambiguity in R-second sampling for GLONASS, with receiver's and satellite's
% differential code biases (RDCB and SDCB) in nanoseconds, elevation angle EL
% in degrees, for site SITE in 4-char, year Y, day-of-year D and hour H, and
% outputs the ambiguity AMB in TECU with the ambiguity-free slant TEC (STEC) in
% TECU for epoches EP in fractional hours.
%
% [EP,STEC,AMB] = AMBIGUITY('SITE',Y,D,[H1 H2],EL,RDCB,SDCB,F,R) specifies the
% period in hours [H1 H2].
%
% F = 1 indicates daily solution; otherwise, hourly solution.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2013/03/21
%   VERSION	:	3.4 2016/12/10

error(nargchk(9,9,nargin));

% Parameters
f1o = 1.602d9; % Base frequency of L1 in Hz
f1i = 0.5625d6; % Frequency increment of L1 in Hz
f2o = 1.246d9; % Base frequency of L2 in Hz
f2i = 0.4375d6; % Frequency increment of L2 in Hz
clight = 2.99792458d8; % Light speed in m/s
ns2m = clight*1e-9; % Conversion factor from nanoseconds to meters
ionfac = 40.3082d16; % Ionospheric factor in m/TECU/s^2
nsat = 26; % Number of GLONASS satellites
re = 6371; % Mean radius of the earth in km
h1 = 650; % Maximum height of ionosphere in km
h2 = 250; % Minimum height of ionosphere in km
rad = 180. /pi;
nepoch = 3600/rate; % Number of epochs in an hour
%minstec = -100; % Minimum tolerance of slant TEC in TECU
%maxstec = 500; % Maximum tolerance of slant TEC in TECU
GNSS = 'R'; % GNSS type

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
l1 = []; l2 = l1; c1 = l1; p2 = l1; epo = l1; % NOTICE: No more use P1

% Load data
obsmatdir =['OBS/glonass/' Y '.' D '/'];
if dflag==1 % Daily
	INFILES = [Site D '0.mat']; % SSSSDDD0.mat
elseif dflag==0 % Hourly
	for i=hr(1):hr(2)-1
		TMPFILE = [Site D char(97+i) '.mat']; % SSSSDDDA.mat
		if exist([obsmatdir TMPFILE],'file')
			INFILES(i-hr(1)+1,:) = TMPFILE; % Collect valid file names
		end
	end
end
for i=1:size(INFILES,1)
	load([obsmatdir INFILES(i,:)],'rinex_version');
	rv = floor(rinex_version);
	switch rv
	case 2
		load([obsmatdir INFILES(i,:)],'C1','P2','L1','L2','year','month','day',...
			'hour','minute','second','event_flag');
		l1 = [l1;L1(~event_flag,:)];
		l2 = [l2;L2(~event_flag,:)];
		c1 = [c1;C1(~event_flag,:)];
		p2 = [p2;P2(~event_flag,:)];

		% Make sure C1 is not all-zero array
		if ~any(any(c1))
			error('All zeros in C1.');
		end

		% Convert 2-digit year into 4-digit year (valid range 1980-2079)
		year = double(year(~event_flag));
		if year>=80
			year = year+1900;
		else
			year = year+2000;
		end
		epo = [epo;datenum(year,double(month(~event_flag)),double(day(~event_flag)),double(hour(~event_flag)),double(minute(~event_flag)),double(second(~event_flag)))]; % Epoches in frational days
		clear C1 L[12] P[12] year month day hour minute second event_flag;
	case 3
		load([obsmatdir INFILES(i,:)],'C1C','L1C','C2P','L2P','year','month',...
			'day','hour','minute','second','event_flag');
		l1 = [l1;L1C(~event_flag,:)];
		l2 = [l2;L2P(~event_flag,:)];
		c1 = [c1;C1C(~event_flag,:)];
		p2 = [p2;C2P(~event_flag,:)];
		epo = [epo;datenum(double(year(~event_flag)),double(month(~event_flag)),double(day(~event_flag)),double(hour(~event_flag)),double(minute(~event_flag)),double(second(~event_flag)))]; % Epoches in fractional days
                clear C1C C2P L1C L2P year month day hour minute second event_flag;
	otherwise
		 error('The RINEX version of the input file should be 2.00 or later.');
	end
end

% Allocate space for observational arrays
ep = (hr(1):rate/3600:hr(2)-rate/3600)'; % New epochs in hours
nep = size(ep,1);
L1 = zeros(nep,nsat); L2 = L1; C1 = L1; P2 = L1;

% Pre-process the observation data
epo = (epo-d0)*24; % Original epoches in hours
n = round((epo-hr(1))*nepoch)+1; % Indices for 1-s or 30-s sampling
L1(n,:) = full(l1); L1(~L1) = NaN;
L2(n,:) = full(l2); L2(~L2) = NaN;
C1(n,:) = full(c1); C1(~C1) = NaN;
P2(n,:) = full(p2); P2(~P2) = NaN;
clear c1 l[12] p[12];

% Remove the data of invalid elevation angles
if any(any(el))
	el = full(el);
	l = ~el;
	L1(l) = NaN; L2(l) = NaN; C1(l) = NaN; P2(l) = NaN; el(l) = NaN;
else
	error('All zeros in variable el.');
end

% Fetch frequency table
switch rv
case 2
	INFILE = ['brdc' num2str(doy,'%03u') 'g.mat'];
	load(['EPM/' INFILE],'sid','freqn');
case 3
	load(INFILES(end,:),'sid','freqn');
end
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
lambda1 = clight./f1; % L1 wavelength in m
lambda2 = clight./f2; % L2 wavelength in m
nr = size(P2,1);
ones_p2 = ones(nr,1);
m2tecu = ones_p2*m2tecu; % Broadcast vector (1 x 26) into array (n x 26)
lambda1 = ones_p2*lambda1;
lambda2 = ones_p2*lambda2;

% Estimate code and phase TEC respectively
stecp = (P2-C1+(ones_p2*glodcb+recdcb)*ns2m).*m2tecu; % Slant TEC based on code observables
stecl = (lambda1.*L1-lambda2.*L2).*m2tecu; % Slant TEC based on phase observables

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
