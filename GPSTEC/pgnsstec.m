%PGNSSTEC   GNSS-TEC parallel processor for Octave.

%   AUTHOR	:	Min-Yang Chou, Ho-Fang Tsai
%   SINCE	:	2016/08/23
%   VERSION	:	0.5 2016/12/15

% Check Octave environment and load parallel package if not loaded
if exist('OCTAVE_VERSION','builtin')
	if ~[pkg('list','parallel'){1}.loaded]
		pkg load parallel;
	end
else
	error('Sorry, Octave only!');
end

% Parameters (change it if necessary)
ff = @bdstec; % GNSS-TEC data processor, which can be gpstec, glotec, galtec or bdstec
Sites=['anmg';'brun';'cmum';'cpnm';'cuut';'dae2';'dltv';'eusm';'hksl';'hkws';'jfng';'jnav';'ncku';'pngm';'sin1'];
%Sites =['banc';'chen';'chia';'henc';'hsin';'hual';'lany';'pang';'suao';'sun1']; % List of GNSS sites
yr = 2016; % Year
doy = 190:192; % Days of year
hrng =[0 24]; % Hour range
sampl = 30; % Sampling period in seconds
minel = 20; % Minimum elevation angle of satellites in degrees
rotate = 1; % Rotate by 'site' or 'doy' (1: date rotating with site parallelization; 2: site rotating with date parallelization)

% Parallel processing
switch rotate
case 1 
	% Run the first case to generate orbital MAT-file (prevent from error)
	n = size(doy,2); % Number of DOY
	Sites1 = mat2cell(reshape(Sites(1,:)',1,4),1,4);
	hrng1  = mat2cell(hrng,1,2);
	for i=1:n
		parcellfun(1,ff,Sites1,{yr},{doy(i)},hrng1,{sampl},{minel});
	end

	% Run the rest
	[nr,nc] = size(Sites(2:end,:)); % Number of sites and number of char. of site name
	prnoc = nr; % Number of processes
	Sites1 = mat2cell(reshape(Sites(2:end,:)',1,nr*nc),1,ones(1,nr)*nc);
	year1  = mat2cell(repmat(yr   ,1,nr),1,ones(1,nr));
	hrng1  = mat2cell(repmat(hrng ,1,nr),1,ones(1,nr)*2);
	sampl1 = mat2cell(repmat(sampl,1,nr),1,ones(1,nr));
	minel1 = mat2cell(repmat(minel,1,nr),1,ones(1,nr));
	for i=1:n
		doy1 = mat2cell(repmat(doy(i),1,nr),1,ones(1,nr));
		parcellfun(prnoc,ff,Sites1,year1,doy1,hrng1,sampl1,minel1);
	end
case 2
	[nr,nc] = size(Sites);
	n = size(doy,2); % Number of DOY
	prnoc = n; % Number of processes
	year1  = mat2cell(repmat(yr   ,1,n),1,ones(1,n));
	doy1   = num2cell(doy);
	hrng1  = mat2cell(repmat(hrng ,1,n),1,ones(1,n)*2);
	sampl1 = mat2cell(repmat(sampl,1,n),1,ones(1,n));
	minel1 = mat2cell(repmat(minel,1,n),1,ones(1,n));
	for i=1:nr
		Sites1 = mat2cell(repmat(Sites(i,:),1,n),1,ones(1,n)*nc);
		parcellfun(prnoc,ff,Sites1,year1,doy1,hrng1,sampl1,minel1);
	end
otherwise
	error('The rotate options can only be either 1 (date rotate) or 2 (site rotate).');
end

