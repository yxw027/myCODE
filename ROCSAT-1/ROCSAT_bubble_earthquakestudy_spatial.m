%%
clear
clc

addpath('E:\AOGS_study')
addpath('E:\Matlabtoolbox')
a = load('eqlist_ROCSAT_20210712_Global3');
earthquake = a.data2;

input = 'E:\Satellite\ROCSAT-1\analyzed\';
output = 'E:\AOGS_study\ROCSAT-1\spatial_analysis\';

for i = 1: 11
    eval(['u' num2str(i) '= find(earthquake(:,10) >= ' num2str(7 + 0.1*(i-1)) ');'])
    eval(['data' num2str(10*(7 + 0.1*(i-1))) ' = earthquake(u' num2str(i) ',:);'])
end

lght = [size(data70,1), size(data71,1), size(data72,1), size(data73,1), size(data74,1),...
    size(data75,1), size(data76,1), size(data77,1), size(data78,1), size(data79,1),...
    size(data80,1)];

datastr = {'data70', 'data71','data72', 'data73', 'data74', 'data75', 'data76',...
    'data77', 'data78', 'data79', 'data80'};

for i = 1: 11
    A = datastr{i};
    eval(['F1_epdata = cell(size(' A ',1), 41);'])
    for j = 1: lght(i)
        disp(j)
        eval(['date = ' A '(j, 1: 3);'])
        
        eval(['eplat = ' A '(j,7);']);
        eval(['eplon = ' A '(j,8);']);
        eval(['epM = ' A '(j,10);']);
        
        doy = day(datetime(date), 'dayofyear');
        eval(['yr = ' A '(j,1);'])
        
        count = 0;
        for k = doy-20: doy+20
            count = count + 1;
            disp(['event_' num2str(j, '%03d') '_day_' num2str(count, '%02d')])
            
            if k < 0 && mod(yr-1,4) ~=0
                fn = dir([input num2str(yr-1) num2str(k+365, '%03d') '.mat.mat']);
            elseif k < 0 && mod(yr-1,4) ==0
                fn = dir([input num2str(yr-1) num2str(k+366, '%03d') '.mat.mat']);
            elseif k > 365 && mod(yr,4) ~=0
                fn = dir([input num2str(yr+1) num2str(k-365, '%03d') '.mat.mat']);
            elseif k > 366 && mod(yr,4) == 0
                fn = dir([input num2str(yr+1) num2str(k-366, '%03d') '.mat.mat']);
            else
                fn = dir([input num2str(yr) num2str(k, '%03d') '.mat.mat']);
            end
            
            if isempty(fn); continue; end
            load([input fn.name])
            GLAT = data.GLAT;
            GLON = data.GLON;
            LogN = data.LogN;
            HHT = data.htpower;
            sigma = logscale.sigma;
            
            R_strain = 10.^(0.43*7);
            GLON(GLON>180) = GLON(GLON>180) - 360;
            
            MAP = zeros(5, 5);

            loncount = 0;
            for LON = eplon-40: 20: eplon+40
                loncount = loncount + 1;
                latcount = 0;
                for LAT = eplat-40: 20: eplat+40
                    disp([num2str(LON) ' ' num2str(LAT)])
                    latcount = latcount + 1;
                    
                    I = find(GLAT >= LAT - 10 & GLAT <= LAT + 10 &...
                        GLON >= LON - 10 & GLON <= LON + 10);
                    
                    if isempty(I);continue;end
                    obsvtn = sqrt(HHT(I))/2;
                    u = find(obsvtn > 4.2);
                    MAP(latcount, loncount) = numel(u);
                end
            end
            
            F1_epdata{j, count} = MAP;
        end
    end
    save([output A ],  A , 'F1_epdata')
end

%%
clear
close all
clc

input = 'E:\AOGS_study\ROCSAT-1\spatial_analysis\';
fn = dir([input '*mat']);

load Mycolormap

figure
doy = -20: 1: 20;
count = 0;
for i = size(fn, 1):-1:8 
    count = count + 1;
    load([input fn(i).name])
    epM = fn(i).name;
    epM = [epM(5) '.' epM(6)];
    
    before = F1_epdata(:, 1: 2);
    after = F1_epdata(:, 40: 41);
    
    b18 = zeros(5, 5, size(F1_epdata,1));
    for j = 1: size(before, 1)
        if isempty(before{j, 2}); continue; end
        b18(:,:,j) = before{j, 2};
    end
    b18 = mean(b18,3);
    
    b19 = zeros(5, 5, size(F1_epdata,1));
    for j = 1: size(before, 1)
        if isempty(before{j, 1}); continue; end
        b19(:,:,j) = before{j, 1};
    end
    b19 = mean(b19,3);
    
    
    a18 = zeros(5, 5, size(F1_epdata,1));
    for j = 1: size(after, 1)
        if isempty(after{j, 1}); continue; end
        a18(:,:,j) = after{j, 1};
    end
    a18 = mean(a18,3);
    
    a19 = zeros(5, 5, size(F1_epdata,1));
    for j = 1: size(after, 1)
        if isempty(after{j, 2}); continue; end
        a19(:,:,j) = after{j, 2};
    end
    a19 = mean(a19,3);
    
    
    MAP = cell(size(after,1), 2);
    for k = 1: 2
        for j = 1: size(after, 1)
            a = after{j, k};
            b = before{j, k};
            logic = b > a;
            MAP{j, k} = logic*1; 
        end
    end
    
    PEIA = zeros(5, 5);
    con = 0;
    for k = 1: 2
        for j = 1: size(MAP, 1)
            if isempty(MAP{j, k}); continue; end
            PEIA = PEIA + MAP{j, k};
            con = con+1;
        end
    end
    numbr = size(MAP, 1);
    PEIA = PEIA/con;
    

    p = PEIA;
    q = 1-p;
    sqrtn = sqrt(numbr);
    denominator = p.*q;
    numerator = sqrtn.*abs((PEIA-0.5));
    Z = numerator./denominator;
    
    index = nan(5, 5);
    index2 = nan(5, 5);
    for j = 1: 5
        
        extP = p(j, :);
        extZ = Z(j, :);
        
        u1 = find(extZ > 2.58 & extP > 0.5);
        u2 = find(isinf(extZ) & extP > 0.5);
        u = union(u1, u2);
        
        u3 = find(extZ > 1.96 & extP > 0.5);
        
        index(j, u) = 1;
        index2(j, u3) = 1;
        
    end
    
    u4 = find(~isnan(index));
    index2(u4) = nan;
    
    
    axes('Position',[0.2 0.75-0.23*(count-1) 0.2 0.20])
    imagesc(-40:20:40, -40:20:40, b18+b19)
    set(gca, 'Ydir', 'normal')
    colormap(flipud(gray))
    caxis([0 3e3])
    axis equal
    xtickangle(0);
    axis([-50 50 -50 50])
    if i ~= 8
        set(gca, 'xtick', [])
        set(gca, 'ytick', -50:20:50)
    else
        set(gca, 'xtick', -50:20:50, 'ytick', -50:20:50)
        xlabel('Days to Earthquake')
        ylabel('Latitude to Epicenter')
        xlabel('Longitude to Epicenter')
    end
    set(gca, 'fontsize', 13)
    if count == 1
        title('Before Earthquake')
    end
    text(-150, 0,[num2str(numbr) '  M > ' epM], 'fontsize', 20)
    
    
    axes('Position',[0.35 0.75-0.23*(count-1) 0.2 0.2])
    imagesc(-40:20:40, -40:20:40, a18+a19)
    set(gca, 'Ydir', 'normal')
    colormap(flipud(gray))
    caxis([0 3e3])
    axis equal
    xtickangle(0);
    axis([-50 50 -50 50])
    if i ~= 8
        set(gca, 'xtick', [])
        set(gca, 'ytick', -50:20:50)
    else
        cb = colorbar('position', [0.52 0.1 0.005 0.75]);
        set(gca, 'xtick', -50:20:50, 'ytick', -50:20:50)
        ylabel(cb, 'Irregularity Events')
    end
    set(gca, 'fontsize', 13)
    if count == 1
        title('After Earthquake')        
    end
    
    
    axes('Position',[0.55 0.75-0.23*(count-1) 0.2 0.2])
    imagesc(-40:20:40, -40:20:40, PEIA)
    set(gca, 'Ydir', 'normal')
    ax = gca;
    colormap(ax, mycolor)   
    caxis([0 1])
    xtickangle(0);
    axis equal
    axis([-50 50 -50 50])
    if i ~= 8
        set(gca, 'xtick', [])
        set(gca, 'ytick', -50:20:50)
    else
        cb = colorbar('position', [0.72 0.1 0.005 0.75]);
        set(gca, 'xtick', -50:20:50, 'ytick', -50:20:50)
        ylabel(cb, 'Proportion')
    end
    set(gca, 'fontsize', 13)
    hold on
    
    
    M = -40: 20: 40;
    count2 = 0;
    for lat = -40: 20: 40
        count2 = count2 + 1;
        P1 = plot(M, lat.*index(count2,:), 'ko', 'markersize', 6, 'markerface', 'k');
        P2 = plot(M, lat.*index2(count2,:), 'ko', 'markersize', 6);
    end
    legend([P1(1), P2(1)], {'P < 0.01', 'P < 0.05'},...
        'Position',[0.84 0.027 0.067 0.069])
    
end

