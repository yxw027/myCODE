function [vo,s] = smooth(v,n,METHOD)
%SMOOTH   Smooth data with running mean.
% [VO,S] = SMOOTH(VI,N) smooths the input vector or array VI with N-point
% running mean along the columns, and returns smoothed vector or array VO and
% its standard deviation.
%
% [VO,S] = SMOOTH(VI,N,'MEDIAN') processes with running median instead of
% running mean, and returns smoothed vector or array and its median absolute
% deivation.

%   AUTHOR	:	Ho-Fang Tsai
%   SINCE	:	2012/04/13
%   VERSION	:	1.3 2016/12/24

error(nargchk(2,3,nargin));

% Convert the input matrix into full matrix if sparsed
sp = issparse(v);
if sp
	v = full(v);
end

% Convert any non-numeric elements into 0.
l = ~isfinite(v);
aal = any(any(l)); % Retrun 1 if any elements are not number
if aal
	v(l) = 0;
end

% Define the critical number
[nr,nc] = size(v);
cn2 = ceil(n/2); % The smallest integer not less than n/2 (2 if n=3 or n=4)
fn2 = floor(n/2); % The largest integer not greater than n/2 (2 if n=4 or n=5)
c8 = round(n*0.8); % Critical number of 80% the given points

% Define the offset
offset = 0;
if cn2==fn2
	offset = offset+1;
end

% Initialize variables
vi = zeros(nr,nc);
s = vi;

% Start N-point running mean/median
if ~exist('METHOD','var')
	for i=1:nc
		for j=cn2+offset:nr-fn2
			v1 = v(j-fn2:j+fn2,i);
			f = find(v1);
			if size(f,1)>c8
				vi(j,i) = mean(v1(f));
				s(j,i) = std(v1(f));
			end
		end
	end
elseif exist('METHOD','var') & lower(METHOD)=='median'
	coef = 1.4826; % Coefficient of median absolute deviation
	for i=1:nc
		for j=cn2+offset:nr-fn2
			v1 = v(j-fn2:j+fn2,i);
			f = find(v1);
			if size(f,1)>c8
				vi(j,i) = median(v1(f));
				s(j,i) = median(abs(v1(f)-vi(j,i)))*coef; % Median absolute deviation (MAD; cf. Leys et al., Detecting outliers: Do not use standard deviation around the mean, use absolute deviation around the median, J. Exp. Soc. Psychol., 49(4), 764-766, 2013.)
			end
		end
	end
else
	warning('METHOD can only be ''median''.');
	return;
end
vo = vi;

% Sparse the output matrix if the input matrix is sparsed; Convert zero
% elements into NaN's if the input matrix contains any non-numeric elements
if sp
	vo = sparse(vo);
elseif aal
	vo(~isfinite(vo)) = NaN;
end

