function crx2rnx(Site,yr,doy,inpath)
%CRX2RNX  Convert daily compact RINEX OBS files into RINEX files.
% CRX2RNX('SITE',Y,D) converts one-day worth of the zipped or unzipped compact
% RINEX OBS files in 30-second sampling into the RINEX files in the same
% sampling rate for site SITE in 4-char, year Y, and day-of-year D.
%
% The source file should be named in SSSSDDD0.YYd, SSSSDDD0.YYd.Z or
% SSSSDDD0.YYd.gz, where DDD is 3-digit day of year, YY is 2-digit year.  The
% output name is in format SSSSDDD0.YYo.
%
% CRX2RNX('SITE',Y,[D1 D2]) specifies the day period [D1 D2].
%
% A site list for multiple sites is also acceptable.  For example,
%
%   SITE = ['aaaa';'bbbb'];
%   CRX2RNX(SITE,2016,1,2)
%
% converts files for both sites.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2015/08/04
%   VERSION	:	0.3 2016/12/14

error(nargchk(4,4,nargin));

% Parameters
Y = num2str(yr,'%04u'); Y2 = Y(3:4);

% Start daily looping
for d=doy(1):doy(end)
	D = num2str(d,'%03u'); % Day of year in 3-char
	for s=1:size(Site,1)

		% Check if the file name exists
		INFILE = [ Site(s,:) D '0.' Y2 'd']; % ssssddd0.yyd
		OUTFILE = [INFILE(1:11) 'o']; % ssssddd0.yyo
		if exist([inpath INFILE],'file')
			system(['crx2rnx -f ' inpath INFILE]);
			disp(['FILE CREATION: ' OUTFILE]);
		elseif exist([inpath INFILE '.Z'],'file') % ssssddd0.yyd.Z
			system(['zcat ' inpath INFILE '.Z | crx2rnx ->' inpath OUTFILE]);
			disp(['FILE CREATION: ' OUTFILE]);
		elseif exist([inpath INFILE '.gz'],'file') % ssssddd0.yyd.gz
			system(['zcat ' inpath INFILE '.gz | crx2rnx ->'  inpath OUTFILE]);
			disp(['FILE CREATION: ' OUTFILE]);
		else
			error([INFILE '[.Z|.gz] not found.']);
		end
	end
end

