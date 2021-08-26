function cof_P=Get_coef(b,s,order)
%GET_COEF   Estimate coefficients of spherical harmonic function.
% COEF = GET_COEF(B,S,ORDER) estimates the ionosphere model coefficients COEF
% of the spherical harmonic basis function, with the given geocentric latitude
% (B) and sun-fixed longitude (S) of the ionosphere pierce point and the order
% of the function.

%   AUTHOR	:	Min-Yang Chou
%   SINCE	:	2015/12/15
%   VERSION	:	0.2 2015/12/18
 
cof_P = zeros(length(b),(order+1)^2);
ms = zeros(length(s),4);
for i=1:length(s)
	ms(i,:)=linspace(s(i),order*s(i),order);
end
i = 1;
x = sin(b);
for n=0:order
	P = legendre(n,x)';
	for m=0:n
		if m==0
			cof_P(:,i) = P(:,m+1).*normP(n,m); % a_n0
		else
			cof_P(:,i) = P(:,m+1).*normP(n,m).*cos(ms(:,m)); % a_nm
			i = i+1;
			cof_P(:,i) = P(:,m+1).*normP(n,m).*sin(ms(:,m)); % b_nm
		end
		i = i+1;
	end
end
