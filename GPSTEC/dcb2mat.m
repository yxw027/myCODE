function dcb2mat(INFILE)
%DCB2MAT  Convert DCB file into MAT-file.
% DCB2MAT('INFILE') converts the DCB file INFILE into a MAT-file.  For example,
%
%  DCB2MAT('P1P21103_ALL.DCB')
%
% converts P1P21103_ALL.DCB into P1P21103_ALL.MAT.
%
%  DCB2MAT('DCB/STATION.DCB')
%
% converts STATION.DCB in the directory DCB into STATION.MAT in the current
% directory (not in DCB/).

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/13
%   VERSION	:	1.6 2015/02/05

error(nargchk(1,1,nargin));

% Load DCB file
[fid, Msg] = fopen(INFILE,'r');
if fid==-1
	error(Msg);
end

% Skip the header of the file
for i=1:7
	fgetl(fid);
end

% Read in the body of the file
i = 1;
while ~feof(fid)
	S = fgetl(fid);
	if size(S,2)<31 % Too few words read
		break;
	end
	PRN(i,1:3) = S(1:3); % GPS PRN or GLONASS slot number;  NOTICE: it would be only 'G' or 'R' if DCBs for stations are read
	StationName(i,1:16) = S(7:22); % Station Name
	dcb(i) = str2double(S(27:min(size(S,2),35))); % DCB in ns
	if size(S,2)>=43
		rms(i) = str2double(S(39:min(size(S,2),47))); % RMS in ns
	else
		rms(i) = NaN;
	end
	i = i+1;
end
fclose(fid);

% Write to file
if strncmp(INFILE,'DCB/',4)
	OUTFILE = [INFILE(5:end-4) '.mat']; % Trim out 'DCB/' and replace the suffix
else
	OUTFILE = [INFILE(1:end-4) '.mat']; % Replace the suffix
end
save(['DCB/' OUTFILE],'PRN','StationName','dcb','rms');
disp(['FILE CREATION: ' OUTFILE]);

