clear
clc


addpath('E:\Matlabtoolbox');
input = 'E:\Satellite\ROCSAT-1\analyzed\';

Year = 2000;
Month = 5;
Date = 4;
DayofYear = day(datetime(Year, Month, Date), 'dayofyear');
suspicious = 0;

epicenter = [-1.12, 123.57];
TsunamiTime = 4 + 21/60 + 16/3600;  % UTC

%--------------------------------------------------------------------------

fn = dir([input num2str(Year) num2str(DayofYear, '%03d') '*']);
load(fn.name)
data.GLON(data.GLON>180) = data.GLON(data.GLON>180) - 360;

index = find(data.Time > TsunamiTime*3600 &...
    data.Time < TsunamiTime*3600+7200);

figure
plotmap

scatter(data.GLON(index), data.GLAT(index), 50, data.LogN(index),'.')

plot(epicenter(2), epicenter(1), 'kp', 'markersize', 20,...
    'markerface', 'k')

xlim([110 130])
ylim([-5 5])
set(gcf, 'color', 'w')
set(gca, 'xtick', -180:10:180, 'ytick', -40:5:40)
caxis([5 6.5])
colormap hsv
cb = colorbar;
ylabel(cb, 'Log(N)')

LON = data.GLON(index);
LAT = data.GLAT(index);
t = data.Time(index);
N = data.LogN(index);



%%
ttt = suspicious;

figure
plot(t(ttt-500:ttt+500), 10.^N(ttt-500:ttt+500))
hold on
line([t(ttt) t(ttt)],[1e6 2e6], 'color', [1 0 0])
set(gca, 'fontsize', 20)
set(gcf, 'color', 'w')
xlabel('Second of Day')
ylabel('Ion Denstiy (#/cm^3)')


r = sqrt(((LAT(ttt)-epicenter(1))*110)^2 + ((LON(ttt)-epicenter(2))*100)^2);
tt = t(ttt)-TsunamiTime*3600-15*60;
v = 1000*r/tt

hold on

plot(LON(2568), LAT(2568), 'bp')


