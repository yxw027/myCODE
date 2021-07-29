clear
clc
input = 'E:\Satellite\ROCSAT-1\analyzed\';
fn = dir([input '*mat']);
Ni = cell(size(fn,1),1);
Ai = cell(size(fn,1),1);
sigma = cell(size(fn,1),1);
for i = 1: size(fn,1)
load([input fn(i).name])
ut = data.Time/3600;
lon = data.GLON;
lon(lon>180) = lon(lon>180) - 360;
LT = ut + lon/15;
LT(LT>24) = LT(LT>24) - 24;
LT(LT<0) = LT(LT<0) + 24;
u1 = find(LT>=20);
u2 = find(LT<2);
u = union(u1,u2);
Ni{i} = data.LogN(u);
Ai{i} = sqrt(data.htpower(u));
sigma{i} = logscale.sigma(u);
disp(i)
end
Ni = cell2mat(Ni);
Ai = cell2mat(Ai);
sigma = cell2mat(sigma);
plot(sigma, log(Ai), 'k.', 'markersize', 0.5)
plot(log(sigma), log(Ai), 'k.', 'markersize', 0.5)
plot(log(sigma*100), log(Ai), 'k.', 'markersize', 0.5)
u = find(sigma*100>0.3);
plot(log(sigma(u)*100), log(Ai(u)), 'k.', 'markersize', 0.5)
plot(log(sigma(u)*100), log(Ni(u)), 'k.', 'markersize', 0.5)
plot(log(sigma(u)*100), Ni(u), 'k.', 'markersize', 0.5)
plot(log(sigma(u)*100), log(Ai(u)), 'k.', 'markersize', 0.5)
plot(Ni, log(Ai), 'k.', 'markersize', 0.5)
loglog(10.^Ni, Ai, 'k.', 'markersize', 0.5)
hold on
line([10^0, 10^7], [10^3.5, 10^3.5], 'color', 'r', 'linewidth', 3)
line([10^0, 10^7], [10^3.6, 10^3.6], 'color', 'r', 'linewidth', 3)
10^3.6
10^3.5
10^3.7
plot(10^3.5, 'r')
line([10^0, 10^7], [10^3.7, 10^3.7], 'color', 'r', 'linewidth', 3)
loglog(10.^Ni, Ai, 'k.', 'markersize', 0.5)
hold on
plot(10^3.5, 'r')
line([10^0, 10^7], [10^3.7, 10^3.7], 'color', 'r', 'linewidth', 3)
plot(log(sigma(u)*100), log(Ai(u)), 'k.', 'markersize', 0.5)
loglog(sigma(u)*100, Ai(u), 'k.', 'markersize', 0.5)