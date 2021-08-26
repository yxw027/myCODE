function dcb = readbia(Type,yr,doy)
%READBIA   Read BIA file and retrieve DCB values for Galileo.
% DCB = READBIA('TYPE',Y,D) reads the BIA file (DLR0MGXFIN_yyyyddd0000_???_07D
% _DCB.BSX) in SINEX format with the given type, say, C1X-C5X, 4-digit year and
% 1- to 3-digit days of the year, and returns the DCB values in nanoseoncds for
% Galileo.  Zeros are returned if values are not available.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2015/09/08
%   VERSION	:	1.1 2017/02/16

error(nargchk(1,3,nargin));

% Parameters
nsat = 30; % Number of satellites
offset = 0; % Offset in columns
dcb = zeros(1,nsat); % Preset DCB vector

% Determine the input file name
switch yr
case 2013, INFILE = 'DCB/DLR0MGXFIN_20130010000_01Y_07D_DCB.BSX';
case 2014, INFILE = 'DCB/DLR0MGXFIN_20140010000_01Y_07D_DCB.BSX';
case 2015, INFILE = 'DCB/DLR0MGXFIN_20150010000_01Y_07D_DCB.BSX';
case 2016,
	if doy<92
		INFILE = 'DCB/DLR0MGXFIN_20160010000_03L_07D_DCB.BSX';
	elseif doy<183
		INFILE = 'DCB/DLR0MGXFIN_20160920000_03L_07D_DCB.BSX';
		offset = -5; % Offset in columns
	elseif doy<275
		INFILE = 'DCB/DLR0MGXFIN_20161830000_03L_07D_DCB.BSX';
		offset = -5;
	else
		INFILE = 'DCB/DLR0MGXFIN_20162750000_03L_07D_DCB.BSX';
		offset = -5;
	end
case 2017, INFILE = 'DCB/DLR0MGXFIN_20162750000_03L_07D_DCB.BSX'; doy = 365;
	offset = -5; % Offset in columns
otherwise, return;
end

% Read in the specified file
[fid, Msg]  = fopen(INFILE,'r');
if fid==-1
	error(Msg);
end

% Read in the file
S = fgetl(fid);
if size(S,2)>=5 && ~strcmp(S(1:5),'%=BIA')
	error('Ths input file is not in BIA format.');
end

% Skip the remarks and unwanted sections
S = fgetl(fid);
while(S(1)=='*')
	S = fgetl(fid);
end
while(size(S,2)~=14 || ~strcmp(S(1:14),'+BIAS/SOLUTION'))
	S = fgetl(fid);
end

% Read in the BIAS section
while(~strcmp(S(1:14),'-BIAS/SOLUTION'))
	S = fgetl(fid);
	if size(S,2)>=92+offset & (strcmp(S(2:4),'DCB') | strcmp(S(2:4),'DSB')) & S(7)=='E' & strcmp(S((31:33)+offset),Type(1:3)) & strcmp(S((36:38)+offset),Type(5:7))
		bst = str2double(S((44:46)+offset)); % Bias start time
		if bst<=doy
			prn = str2double(S(13:14)); % PRN
			if strcmp(S(13:14),'  ')
				break;
			end
			dcb(prn) = str2double(S((72:92)+offset)); % DCB for the given type
		end
	end
end
fclose(fid);

