function flag = checkvar(INFILE,V)
%CHECKVAR  Check if MAT-file contains specified variable.
% FLAG = CHECKVAR('FILE','VAR') checks if the MAT-file FILE contains variable
% VAR, and returns 1 if true; otherwise, 0.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2016/08/20
%   VERSION	:	0.2 2016/08/28

% Parameters
flag = false; % Preset flag

% Start checking
C = who('-file',INFILE);
for i=1:size(C,1)
	if strcmp(C{i},V)
		flag = true;
		return;
	end
end
