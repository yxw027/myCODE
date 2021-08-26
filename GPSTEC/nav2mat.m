function nav2mat(INFILE)
%NAV2MAT  Convert a RINEX NAV file into a MAT-file for GLONASS.
% NAV2MAT('SSSSDDD0.YYg') converts a GLONASS NAV file, SSSSDDD0.YYg, in RINEX
% format into a MAT-file.  For example, NAV2MAT('brdc0010.12g') converts
% 'brdc0010.12g' into 'brdc001g.mat'.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/10/11
%   VERSION	:	0.9 2016/08/23

error(nargchk(1,1,nargin));

% Parameters
if INFILE(8)~='0' % Daily file index
	error('Only one day worth of the NAV file, i.e, SSSSDDD0.YYg can be processed.');
end
if INFILE(12)~='g' % GNSS type
	error('Only RINEX NAV file for GLONASS, i.e., SSSSDDD0.YYg, can be processed.');
end

% Open file
EPMPATH =['EPM/'];
[fid, Msg] = fopen([EPMPATH INFILE],'r');
if fid==-1
	error(Msg);
end

% Skip the header of the file
while ~feof(fid)
	S = fgetl(fid);
	S(size(S,2)+1:80)=char(32); % Compensate the tail with whitespaces if necessary
	if strcmp(S(61:73),'END OF HEADER')
		break;
	end
end

% Read in the body of the file
i = 1;
while ~feof(fid)
	% Read the first line
	S = fgetl(fid);
	SN(i,:)   = S(1:2); % Satellite (slot) numbers in 2-char
	Year(i,:) = S(4:5); % Year in 2-char
	M(i,:)    = S(7:8); % Month in 2-char
	D(i,:)    = S(10:11); % Day in 2-char
	H(i,:)    = S(13:14); % Hour in 2-char
	Mi(i,:)   = S(16:17); % Minute in 2-char
	Se(i,:)   = S(19:22); % Second in %4.1f
	CB(i,:)   = S(23:41); % SV clock bias in seconds
	RFB(i,:)  = S(42:60); % SV relative frequency bias
	MFT(i,:)  = S(61:79); % Message frame time

	% Read the second line
	S = fgetl(fid);
	X(i,:)    = S(4:22); % Satellite position X in km
	Vx(i,:)   = S(23:41); % Velocity X dot in km/s
	Ax(i,:)   = S(42:60); % X acceleration in km s^(-2)
	Heal(i,:) = S(61:79); % Health (0=OK)

	% Read the third line
	S = fgetl(fid);
	Y(i,:)    = S(4:22); % Satellite position Y in km
	Vy(i,:)   = S(23:41); % Velocity Y dot in km/s
	Ay(i,:)   = S(42:60); % Y acceleration in km s^(-2)
	FN(i,:)   = S(61:79); % Frequency number (-7 ... +13) 

	% Read the fourth line
	S = fgetl(fid);
	Z(i,:)    = S(4:22); % Satellite position Z in km
	Vz(i,:)   = S(23:41); % Velocity Z dot in km/s
	Az(i,:)   = S(42:60); % Z acceleration in km s^(-2)
	AOI(i,:)  = S(61:79); % Age of operation information in days
	i = i+1;
end
fclose(fid);

% Convert char into digital
sid	= int8(str2double(SN));
years   = int8(str2double(Year)); % 2-digit year
months  = int8(str2double(M));
days    = int8(str2double(D));
hours   = int8(str2double(H));
minutes = int8(str2double(Mi));
seconds = str2double(Se);
cbias   = str2num(CB);
rfbias  = str2num(RFB);
t_mf    = str2num(MFT);
x       = str2num(X);
y       = str2num(Y);
z       = str2num(Z);
vx      = str2num(Vx);
vy      = str2num(Vy);
vz      = str2num(Vz);
ax      = str2num(Ax);
ay      = str2num(Ay);
az      = str2num(Az);
health  = logical(str2num(Heal));
freqn   = int8(str2num(FN));
age     = str2num(AOI);

% Write to file
OUTFILE = [EPMPATH INFILE(1:7) INFILE(12) '.mat']; % SSSSDDDg.mat
save(OUTFILE,'sid','years','months','days','hours','minutes','seconds',...
	'cbias','rfbias','t_mf','x','y','z','vx','vy','vz','ax','ay','az',...
	'health','freqn','age');
disp(['FILE CREATION: ' OUTFILE]);
