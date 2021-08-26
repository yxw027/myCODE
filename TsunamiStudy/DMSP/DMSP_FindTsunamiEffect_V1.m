clear
clc
close all
warning('off')

addpath('E:\Matlabtoolbox')
addpath('E:\CODE\TsunamiStudy\DMSP')
load hdcoast
load MagLat
ScreenSize = get(0, 'screensize');
Satellite = 'f15';
cd(['E:\Satellite\DMSP\' Satellite])

Epicenter = [38.6, 142.5];
OnsetTime = [5, 46, 23];
EarthquakeTime = OnsetTime(1) + OnsetTime(2)/60 + OnsetTime(3)/3600;
minUT = EarthquakeTime + 0;
maxUT = EarthquakeTime + 2;

Year = 2011;
count = 0;
ObservationDays = 63: 70;

plotmap()
f15data = cell(8, 3);
%--------------------------------------------------------------------------
for Doy = ObservationDays % 99
    
    CumulativeGeoLat = cell(numel(ObservationDays), 1);
    CumulativeGeoLon = cell(numel(ObservationDays), 1);
    CumulativeGeoAlt = cell(numel(ObservationDays), 1);
    CumulativeTime = cell(numel(ObservationDays), 1);
    CumulativePlasmaDensity = cell(numel(ObservationDays), 1);
    CumulativeHorizontalVelocity = cell(numel(ObservationDays), 1);
    CumulativeVerticalVelocity = cell(numel(ObservationDays), 1);
    CumulativeIonDensity = cell(numel(ObservationDays), 1);
    
    Date = datetime(Year, 1, Doy);
    Month = Date.Month;
    Day = Date.Day;
    
    fn = dir(['*2011' num2str(Month, '%02d') num2str(Day, '%02d') '*.mat']);
    %     if isempty(fn);continue;end
    
    disp([num2str(Year, '%02d') '_' num2str(Doy, '%03d')])
    
    data = load(fn.name);
    count = count + 1;
    Time = InterpalTIME(data.Time);
    
    geoLat = nan(size(data.IonDensity, 1), 1); % whole day location
    geoLon = nan(size(data.IonDensity, 1), 1); % whole day location
    geoAlt = nan(size(data.IonDensity, 1), 1); % whole day location
    
    for j = 1: size(data.geoLat, 1)
        a1 = data.geoLat(j, :);
        geoLat(1+20*(j-1), :) = a1;
        a2 = data.geoLon(j, :);
        geoLon(1+20*(j-1), :) = a2;
        a3 = data.geoAlt(j, :);
        geoAlt(1+20*(j-1), :) = a3;
    end
    
    DataLength = (1: 1: size(data.IonDensity, 1))';
    
    geoLat = interp1(DataLength, geoLat, DataLength, 'PCHIP');
    geoLon = interp1(DataLength, geoLon, DataLength, 'PCHIP');
    geoAlt = interp1(DataLength, geoAlt, DataLength, 'PCHIP');
    
    geoLon(geoLon>180) = geoLon(geoLon>180) - 360;
    
    CumulativeTime{count,1} = datevec(Time);
    CumulativeGeoLat{count,1} = geoLat;
    CumulativeGeoLon{count,1} = geoLon;
    CumulativeGeoAlt{count,1} = geoAlt;
    CumulativePlasmaDensity{count,1} = data.PlasmaDensity;
    CumulativeHorizontalVelocity{count,1} = data.HorizontalVelocity;
    CumulativeVerticalVelocity{count,1} = data.VerticalVelocity;
    CumulativeIonDensity{count,1} = data.IonDensity;
    
    CumulativeTime = cell2mat(CumulativeTime);
    CumulativeGeoLat = cell2mat(CumulativeGeoLat);
    CumulativeGeoLon = cell2mat(CumulativeGeoLon);
    CumulativeGeoAlt = cell2mat(CumulativeGeoAlt);
    CumulativePlasmaDensity = cell2mat(CumulativePlasmaDensity);
    CumulativeHorizontalVelocity = cell2mat(CumulativeHorizontalVelocity);
    CumulativeVerticalVelocity = cell2mat(CumulativeVerticalVelocity);
    CumulativeIonDensity = cell2mat(CumulativeIonDensity);
    
    CumulativePlasmaDensity(CumulativePlasmaDensity < 0) = nan;
    CumulativeHorizontalVelocity(CumulativeHorizontalVelocity < -1e36) = nan;
    CumulativeVerticalVelocity(CumulativeVerticalVelocity < -1e36) = nan;
    
    UTHour = CumulativeTime(:, 4) + CumulativeTime(:, 5)/60 + CumulativeTime(:, 5)/3600;
    LTHour = UTHour + CumulativeGeoLon/15;
    LTHour(LTHour<0) = LTHour(LTHour<0) + 24;
    LTHour(LTHour>24) = LTHour(LTHour>24) - 24;
    
    % distinguish
    [OrbitTime, OrbitLat, OrbitLon, OrbitAlt, OrbitPlasmaDensity,...
        OrbitHorizontal, OrbitVertical] = ExtractandDistinguish(CumulativeTime,...
        CumulativeGeoLat, CumulativeGeoLon, CumulativeGeoAlt,...
        CumulativePlasmaDensity, CumulativeHorizontalVelocity,...
        CumulativeVerticalVelocity, Epicenter(2), Epicenter(1), UTHour, minUT, maxUT);
    
%     o = 1;
%     if Doy == 69 || Doy == 70
%         o = 2;
%     elseif Doy == 67 || Doy == 68
%         o = 3;
%     end
     
o = 1;
if Doy == 67
    o = 3;
end

    scatter(OrbitLon{o}, OrbitLat{o}, 100, OrbitPlasmaDensity{o}, '.')
    set(gcf, 'color', 'w')
    caxis([0 5.5e4])
    axis equal
    axis([120 180 20 60])
    hold on
    plot(Epicenter(2), Epicenter(1), 'kp', 'markersize', 12, 'markerface', 'r')
    
    index = OrbitLat{o} > 30 & OrbitLat{o} < 50;
    
    text(median(OrbitLon{o}(index)), median(OrbitLat{o}(index)), num2str(Doy, '%02d'))
    
    f15data{count, 1} = OrbitLon{o};
    f15data{count, 2} = OrbitLat{o};
    f15data{count, 3} = OrbitPlasmaDensity{o};
    
end
cb = colorbar;
ylabel(cb, 'Plasma Density (#/cm^3)');
title('DMSP/f15 DOY 063 - DOY 070 2011')

%%
figure
subplot(411)
plot(orbit{1,1}(:,1), orbit{1,2})
set(gca, 'fontsize', 15)

subplot(412)
plot(orbit{1,3}(:,1), orbit{1,2})
set(gca, 'fontsize', 15)

subplot(413)
plot(orbit{1,3}(:,1), orbit{1,5})
set(gca, 'fontsize', 15)

subplot(414)
plot(orbit{1,3}(:,1), orbit{1,4})
set(gca, 'fontsize', 15)

set(gcf, 'color', 'w')



% u = find(orbit{3,1}(:,1)>15, )

%% Plasma Density

figure('position', [sz(1) sz(2) 0.6*sz(3) 0.6*sz(4)])
subplot('position', [0.25 0.55 0.56 0.4])


L1 = cell2mat(L);
P1 = cell2mat(P);
VH1 = cell2mat(VH);
VV1 = cell2mat(VV);

c = jet(max(P1));
for con = 1: size(P1, 1)
    try
        disp(num2str(con))
        plot(L1(con, 2), L1(con, 1),  '.', 'color', c(floor(P1(con)),:))
        hold on
    catch
        continue
    end
end
title(['DMSP/' sat ' DOY' num2str(doy, '%03d') ', 2011 1547LT-2400LT'] )
plot(hdcoast(:,1), hdcoast(:,2), 'k')
hold on


set(gca,'xtick', -180:60:180, 'ytick', -90:30:90)
axis([-180 180 -90 90])
colorbar
colormap jet

axes('position', [0.242 0.08 0.455 0.857])
set(gca, 'xtick', 0:0.33:1, 'xticklabel', 5:8, 'ytick', [], 'fontsize', 20)
xlabel('UT')

box on
d = 63;
for o = 1: size(L, 1)
    axes('position', [0.25 0.05+0.1*(o-1) 0.5 0.2])
    disp(num2str(o))
    L2 = L{o, 1}; % Location
    P2 = P{o, 1}; % Plasma-
    S2 = S{o, 1};
    VH2 = VH{o, 1};
    VV2 = VV{o, 1};
    %     D = orbit{o, 3}; % Date
    
    plot(S2, VV2)
    hold on
    %         plot(1:2500, 2e4*ones(1,2500),'r--')
    plot(18001:28800, 0*ones(1,10800),'r--')
    %         plot(1:2500, 6e4*ones(1,2500),'r--')
    text(16000, 0, ['DOY ' num2str(d, '%03d')], 'fontsize', 15)
    if o == 1
        text(29000, 0, [num2str(0) ' (#/cm^3)'], 'fontsize', 15, 'color', 'r')
    end
    d=d+1;
    box off
    set(gca,'Visible','off');
end

set(gcf,'color','w');
xlabel('Number (count)')
ylabel('Density  #/c.c.')

close all

arr = 1 + 42/60+ 39/3600;

distant = sqrt(100*((83-epicenter(2)))^2 + (110*(16-epicenter(1)))^2);

plot([83,epicenter(2)], [16, epicenter(1)])

vs = 1000*distant/((arr - earthqk) * 3600)