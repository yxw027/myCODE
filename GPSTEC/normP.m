function N=normP(n,m)
%NORMP   Normalization function of Legendre polynomials
% N = NORMP(N,M) estimates the normalization function of Legendre polynomials,
% N, with the given order n and degree m of the function.

%   AUTHOR	:	Min-Yang Chou
%   SINCE	:	2015/12/15
%   VERSION	:	0.2 2015/12/18

if m==0
	N = sqrt(factorial(n-m)*(2*n+1)/factorial(n+m));
else
	N = sqrt(factorial(n-m)*(4*n+2)/factorial(n+m));
end
