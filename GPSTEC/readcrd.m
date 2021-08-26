function [x,y,z] = readcrd(Site,Domes)
%READCRD   Read CRD files and retrieve coordinate information.
% [X, Y, Z] = READCRD('SITE') reads the coordinate files (currently, CODE.CRD
% and LOCAL.CRD) and retrieves XYZ coordinate values in km for site SITE (4-
% char, case insensetive).
%
% [X, Y, Z] = READCRD('SITE','DOMES') further assigns the DOMES number (format
% nnnnnMnnnXX or nnnnnSnnnXX and case sensetive).
%
% If the specified site not found in the CRD files, site-specific file in form
% of SSSS_nnnnnMnnnXX.mat, SSSS_nnnnnSnnnXX.mat, or SSSS____________.mat (if no
% Domes number specified) are read.  

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2013/08/01
%   VERSION	:	1.0 2016/09/21

narginchk(1,2);

% Parameters
INFILE{1} = 'LOCAL.CRD';
INFILE{2} = 'CODE.CRD';
m2km = 1d-3; % Conversion factor between m and km
r_min = 6.3d3; % Minimum of the possible length of a receiver vector in km
cols = 6:21; % Column number in the CRD file

% Construct full marker name
if size(Site,2)==4
	Site_uc = upper(Site);
else
	error('SITE must be in 4-char.');
end
if exist('Domes','var')
	Domes = deblank(Domes);
	nDomes = size(Domes,2);
	if nDomes>=9 && nDomes<11
		Domes = [Domes char(32*ones(1,11-nDomes))]; % Compensate the tail with whitespaces
	else
		Domes = char(32*ones(1,11)); % All blanks
	end
else
	Domes = char(32*ones(1,11)); % All blanks
end
Site_fn = [Site_uc ' ' Domes]; % Full marker name in 16-char

% Read receiver vector from the CRD files
for j = 1: size(INFILE,2)
    fid = fopen(INFILE{j}, 'r');
    if fid > 0
        for i=1:6
            fgetl(fid); % Skip the header
        end
        while ~feof(fid) % Read XYZ for the matched site name
            S = fgetl(fid);
            if size(S,2)>=66 && strcmp(S(cols),Site_fn)
                x = str2double(S(22:36))*m2km; % X-coord. in km (ECF)
                y = str2double(S(37:51))*m2km; % Y-coord. in km (ECF)
                z = str2double(S(52:66))*m2km; % Z-coord. in km (ECF)
                rn = norm([x y z]); % Length of the vector
                break;
            end
        end
        fclose(fid);
        if exist('rn','var') % Jump out of the loop once XYZ is found
            break;
        end
    else
        warning(['Cannot find ' INFILE{j} '.']);
    end
end

% Read in the site-specific file if XYZ not found
INFILE{3} = [Site_fn '.mat']; % SSSS nnnnn[MS]nnnXX.mat
INFILE{3}(INFILE{3}==' ') = '_'; % SSSS_nnnnn[MS]nnnXX.mat
if ~exist('rn','var')
	if exist(INFILE{3},'file')
		load(INFILE{3},'x','y','z'); % in km (ECF)
		rn = norm([x y z]);
	else
		error(['Cannot find Site ' Site_fn ' in the CRD files.']);
	end
end

% Error if length of 'r' is too short
if rn<r_min
	error(['The length of the receiver vector is too short (|r| = ' int2str(rn) ' km).']);
end

