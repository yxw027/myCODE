clear
clc

addpath('E:\AOGS_study')
a = load('eqlist_20210708_Global3');
earthquake = a.data2;

input = 'E:\Satellite\ROCSAT-1\analyzed\';
output = 'E:\AOGS_study\ROCSAT-1\MAP\';
year = 1999;
lat = 23.85;
lon = 121.82;

LON = -179.5: 1: 179.5;
LAT = 40: -1: -40;
lonr = 2.5;
latr = 2.5;
period = 30;

for dayrange = 1: 30
    
    eqday = datetime(year, 9, dayrange, 17, 47, 0);
    doy = day(eqday, 'dayofyear');
    
    GLON = cell(15, 1);
    GLAT = cell(15, 1);
    LogN = cell(15, 1);
    HHTpower = cell(15, 1);
    Time = cell(15, 1);
    Temperature = cell(15, 1);
    sigma = cell(15, 1);
    deltaN = cell(15, 1);
    
    count = 0;
    for i = doy-floor((period/2)): doy+floor((period/2))
        disp(i)
        count = count + 1;
        try
            load([input num2str(year) num2str(i, '%03d') '.mat.mat']);
        catch
            continue
        end
        
        GLON{count} = data.GLON;
        GLAT{count} = data.GLAT;
        LogN{count} = data.LogN;
        Time{count} = data.Time;
        HHTpower{count} = data.htpower;
        Temperature{count} = data.Temp;
        
        sigma{count} = logscale.sigma;
        
        deltaN{count} = linearscale.deltaNi(1:86400);
    end
    
    GLON = cell2mat(GLON);
    GLON(GLON>180) = GLON(GLON>180) - 360;
    GLAT = cell2mat(GLAT);
    LogN = cell2mat(LogN);
    Time = cell2mat(Time);
    HHTpower = cell2mat(HHTpower);
    Temperature = cell2mat(Temperature);
    sigma = cell2mat(sigma);
    deltaN = cell2mat(deltaN);
    
    UT = Time/3600;
    LT = UT + GLON/15;
    LT(LT>24) = LT(LT>24) - 24;
    LT(LT<0) = LT(LT<0) + 24;
    
    u1 = find(LT>20);
    u2 = find(LT<2);
    u = union(u1, u2);
    
    Ai = mymap(LON, LAT, lonr, latr, GLON, GLAT, sqrt(HHTpower), 'mean');
    
    Si = mysigma(LON, LAT, lonr, latr, GLON, GLAT, sigma, 0.3);
    
    dNi = mymap(LON, LAT, lonr, latr, GLON, GLAT, deltaN, 'mean');
    
    switch period
        case 15
            name = ['ROCSAT-1_chichi_bubble_+-7_09' num2str(dayrange)];
        case 30
            name = ['ROCSAT-1_chichi_bubble_+-15_09' num2str(dayrange)];
    end
    
    save([output name], 'Ai', 'Si', 'dNi')
end