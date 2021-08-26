function sflag = putintv(INFILE,intv)
%PUTINTV   Insert observation interval into RINEX file.
% FLAG = PUTINTV('FILENAME',INTV) inserts observation interval INTV in seconds
% into the header of the specified RINEX file, and returns 1 if complete; 0 if
% fail.  For example,
%
%   PUTINTV('aaaa1680.11d',30)
%
% inserts 30-s interval into the header of file 'aaaa1680.11d'.  If INTERVAL
% already exists, do nothing and return 0.
%
%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2016/12/15
%   VERSION	:	0.2 2016/12/22

error(nargchk(2,2,nargin));

% Parameters
sflag = false;

% Open files
TMPFILE = [INFILE '.tmp'];
fid1 = fopen(INFILE,'r');
fid2 = fopen(TMPFILE,'w');

while ~feof(fid1)
	S = fgetl(fid1);
	if ~sflag && size(S,2)>=77 && strcmp(S(61:77),'TIME OF FIRST OBS')
		fprintf(fid2,...
'%10.3f                                                  INTERVAL\n',intv);
		sflag = true;
	elseif size(S,2)>=68 && strcmp(S(61:68),'INTERVAL')
		sflag = false;
		warning('INTERVAL already exists.');
		fclose(fid2);
		fclose(fid1);
		delete(TMPFILE);
		return;
	end
	fputs(fid2,[S char(10)]);
end
fclose(fid1);
fclose(fid2);
movefile(TMPFILE,INFILE);

