function rinex_rename(Site,yr,doy)
%RINEX_RENAME  Change new RINEX filename parameters to original format.
% RINEX_RENAME('SITE',Y,D) changes file name from new naming convention
% (SSSS?????_?_yyyydddhhmm_01D_30S.??.crx) to the original (ssssddd0.yyd).
%  For example,
%
%   RINEX_RENAME('ptgg',2016,189)
%
% renames PTGG?????_?_20161890000_01D_30S.??.crx into ptgg1890.16d.
%
% RINEX_RENAME('SITE',Y,[D1 D2]) renames files for the day period [D1 D2].
%
% A site list for multiple sites is also acceptable.  For example,
%
%   SITE = ['aaaa';'bbbb'];
%   RINEX_RENAME(SITE,2016,1,2)
%
% renames files for both sites.
%
% The observation interval in seconds will be inserted into the destination
% file if not found in the header of the file.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2016/12/14
%   VERSION	:	0.2 2016/12/22

error(nargchk(3,3,nargin));

% Parameters
Site_lc = lower(Site);
Site_uc = upper(Site);
Y = num2str(yr,'%04u'); % Year in 4-char
Y2 = Y(3:4);

for d=doy(1):doy(end)
	D = num2str(d,'%03u'); % Day-of-year in 3-char
	for s=1:size(Site,1)
		INFILES = dir([Site_uc(s,:) '?????_?_' Y D '0000_01D_30S_??.crx']); % SSSS?????_?_yyyyddd0000_01D_30S_??.crx
		for i=1:size(INFILES,1)
			INFILE = INFILES(i).name;
			OUTFILE = [Site_lc(s,:) D '0.' Y2 'd']; % ssssddd0.yyd
			[err,Msg] = rename(INFILE,OUTFILE);
			if err
				warning(Msg);
			else
				disp(['FILE RENAMING: ' INFILE ' -> ' OUTFILE '.']);
			end

			% Insert observation interval if not found
			sflag = putintv(OUTFILE,str2double(INFILE(29:30)));
		end
	end
end
