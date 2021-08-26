%%
clear
clc

addpath('E:\AOGS_study')
a = load('eqlist_ROCSAT_20210712_Global3');
earthquake = a.data2;

input = 'E:\Satellite\ROCSAT-1\analyzed\';
output = 'E:\AOGS_study\ROCSAT-1\test\';

for i = 1: 16
    eval(['u' num2str(i) '= find(earthquake(:,10) >= ' num2str(6.5 + 0.1*(i-1)) ');'])
    eval(['data' num2str(10*(6.5 + 0.1*(i-1))) ' = earthquake(u' num2str(i) ',:);'])
end

lght = [size(data65,1), size(data66,1), size(data67,1), size(data68,1), size(data69,1),...
    size(data70,1), size(data71,1), size(data72,1), size(data73,1), size(data74,1),...
    size(data75,1), size(data76,1), size(data77,1), size(data78,1), size(data79,1),...
    size(data80,1)];

datastr = {'data65', 'data66', 'data67', 'data68', 'data69', 'data70', 'data71',...
    'data72', 'data73', 'data74', 'data75', 'data76', 'data77', 'data78'....
    , 'data79', 'data80'};

for i = 1: 16
    A = datastr{i};
    eval(['F1_epdata = cell(size(' A ',1), 90);'])
    for j = 1: lght(i)
        disp(j)
        eval(['date = ' A '(j, 1: 3);'])
        
        eval(['eplat = ' A '(j,7);']);
        eval(['eplon = ' A '(j,8);']);
        eval(['epM = ' A '(j,10);']);
        
        doy = day(datetime(date), 'dayofyear');
        eval(['yr = ' A '(j,1);'])
        
        count = 0;
        for k = doy-45:  doy+45
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
            
            u = find(GLAT >= eplat - R_strain/110 & GLAT <= eplat + R_strain/110 &...
                GLON >= eplon - R_strain/105 & GLON <= eplon + R_strain/105);
            
            F1_epdata{j, count} = [GLAT(u), GLON(u), LogN(u), HHT(u), sigma(u)];
        end
    end
    save([output A ],  A , 'F1_epdata')
end

%%
clear
clc

input = 'E:\AOGS_study\ROCSAT-1\global_event\Global_PEIA_range01_lowlat\';
fn = dir([input '*mat']);

sz = get(0, 'screensize');

datastr = {'data65', 'data66', 'data67', 'data68', 'data69', 'data70', 'data71',...
    'data72', 'data73', 'data74', 'data75', 'data76', 'data77', 'data78'....
    , 'data79', 'data80'};

for i = 13
    load([input fn(i).name])
    
    earthquake_bubble = nan(size(F1_epdata));
    earthquake_number = size(F1_epdata, 1);
    for j = 1: earthquake_number
        for k = 1: 91
            
            a = F1_epdata{j, k};
            if isempty(a); continue; end
            u = find(sqrt(a(:,4))>10^4);
            
            earthquake_bubble(j, k) = numel(u);
            
        end
    end
    earthquake_bubble(isnan(earthquake_bubble)) = 0;
end

[pst, index] = sort(data77(:, 10));
data77 = data77(index, :);


figure('position', sz, 'visible', 'off')
for i = 1: 12
    c = i;
    axes('position', [0.14 0.25+(i-1)*0.06 0.10 0.05])
    event = earthquake_bubble(c,:);
    bar(-45:1:45, event)
    box on
    
    set(gca, 'fontsize', 10)
    if i == 1
        h = gca;
        set(gca, 'Ytick', 0:100:500, 'xtick', -45:10:45)
        % rotateticklabel(h,0)
        xtickangle(0);
        xlabel('Day to Earthquake')
    elseif i == 12
        title('Irregularity Events')
        set(gca, 'Ytick', 0:100:500, 'xtick', -45:10:45)
        set(gca, 'xticklabel', [], 'yticklabel', {})
    else
        set(gca, 'Ytick', 0:100:500, 'xtick', -45:10:45)
        set(gca, 'xticklabel', [], 'yticklabel', {})
    end
    
    axis([-45 45 0 500])
    grid on
    box on
    line([0 0], [0 500], 'color', 'red')
    
    lat = num2str(data77(c,7));
    lon = num2str(data77(c,8));
    M = num2str(data77(c, 10));
    strdate = datestr(data77(c,1:6), 26);
    
    text(-150, 320, strdate, 'fontsize', 10)
    text(-110, 320, ['M = ' M], 'fontsize', 10)
    text(-155, 200, ['epicenter = (' lat ', ' lon ')'],...
        'fontsize', 10)
    
%-------------------------------------------------------------------------%
%                      cumulative bubble events                           %
%-------------------------------------------------------------------------%
    
    before = fliplr(event(1:45));
    after = event(47:91);
    
    cumulative_before = nan(1, 45);
    cumulative_after = nan(1, 45);
    for j = 1: 45
        cumulative_before(j) = sum(before(1:j), 'omitnan');
        cumulative_after(j) = sum(after(1:j), 'omitnan');
    end
    
    axes('position', [0.27 0.25+(i-1)*0.06 0.20 0.05])
    hold on
    bar(0.5:1:44.5, cumulative_before, 0.4, 'c')
    bar(1:1:45, cumulative_after, 0.4, 'r')
    box on
    set(gca, 'fontsize', 10)
    
    if i == 1
        h = gca;
        legend({'Before', 'After'}, 'Position',...
            [0.427 0.18 0.04 0.03])
        set(gca, 'Ytick', 0:1000:5000, 'xtick', 0:5:45)
        xlabel('Day to Earthquake')
        xtickangle(0);
    elseif i == 12
        set(gca, 'Ytick', 0:1000:5000, 'xtick', 0:5:45)
        title('Cumulative Irregularity Events')
        set(gca, 'xticklabel', [], 'yticklabel', {})
    else
        set(gca, 'Ytick', 0:1000:5000, 'xtick', 0:5:45)
        set(gca, 'xticklabel', [], 'yticklabel', {})
    end
    axis([0 45 0 5000])
    grid on
    box on
    
%-------------------------------------------------------------------------%
%                   Pre-Earthequake Irregularity Anomaly                  %
%-------------------------------------------------------------------------%
    
    axes('position', [0.50 0.25+(i-1)*0.06 0.10 0.05])
    
    PELA = cumulative_before > cumulative_after;
    plot(PELA, '*', 'markersize', 3)
    axis([0 45 0.9 1.1])
    
    set(gca, 'fontsize', 10)
    if i == 1
        h = gca;
        set(gca, 'xtick', 0:5:45, 'yticklabel', [])
        xlabel('Day to Earthquake')
        xtickangle(0);
    elseif i == 12
        set(gca, 'xtick', 0:5:45)
        title('Pre-Earthequake Irregularity Anomaly')
        set(gca, 'xticklabel', [], 'yticklabel', {})
    else
        set(gca, 'xtick', 0:5:45)
        set(gca, 'xticklabel', [], 'yticklabel', {})
    end
end

axes('position', [0.65 0.91 0.10 0.05])
PELA = nan(size(data77,1), 45);
for i = 1: size(data77, 1)
    
    event = earthquake_bubble(i,:);
    before = fliplr(event(1:45));
    after = event(47:91);
    
    cumulative_before = nan(1, 45);
    cumulative_after = nan(1, 45);
    for j = 1: 45
        cumulative_before(j) = sum(before(1:j), 'omitnan');
        cumulative_after(j) = sum(after(1:j), 'omitnan');
    end
    logic = cumulative_before > cumulative_after;
    PELA(i, :) = logic*1;
end

sumPELA = sum(PELA, 1);
plot(sumPELA)
ylim([0 12])
hold on
line([0 45],[12/2 12/2],'color','r')
set(gca, 'fontsize', 10)
xlabel('Day to Earthquake')
ylabel('Event Count')
set(gca, 'xtick', 0:5:45)
xtickangle(0);

path_f = 'E:\AOGS_study\ROCSAT-1\';
set(gcf,'paperpositionmode','auto')
print('-dpng','-r300',[path_f 'ROCSAT_temporal.jpg'])
close all

%%
clear
clc

input = 'E:\AOGS_study\ROCSAT-1\global_event\Global_PEIA_area5\';
output = 'E:\AOGS_study\ROCSAT-1\global_event\Global_PEIA_area5\';

fn = dir([input '*mat']);

datastr = {'data65', 'data66', 'data67', 'data68', 'data69', 'data70', 'data71',...
    'data72', 'data73', 'data74', 'data75', 'data76', 'data77', 'data78'....
    , 'data79', 'data80'};

sz = get(0, 'screensize');

count = 0;
for power = 3.5: 0.1: 5
    count = count + 1;
    proportion = nan(size(fn,1), 45);
    event_number = nan(size(fn,1),1);
    for i = 1: size(fn, 1)
        load([input fn(i).name])
        A = datastr{i};
        
        earthquake_bubble = nan(size(F1_epdata));
        earthquake_number = size(F1_epdata, 1);
        for j = 1: earthquake_number
            for k = 1: 91
                a = F1_epdata{j, k};
                if isempty(a); continue; end
                u = find(sqrt(a(:,4)) > 10^power);
                earthquake_bubble(j, k) = numel(u);
            end
        end
        earthquake_bubble(isnan(earthquake_bubble)) = 0;
        
        eval(['PELA = nan(size(' A ',1), 45);'])
        for j = 1: size(earthquake_bubble, 1)
            
            event = earthquake_bubble(j,:);
            before = fliplr(event(1:45));
            after = event(47:91);
            
            cumulative_before = nan(1, 45);
            cumulative_after = nan(1, 45);
            for k = 1: 45
                cumulative_before(k) = sum(before(1:k), 'omitnan');
                cumulative_after(k) = sum(after(1:k), 'omitnan');
            end
            logic = cumulative_before > cumulative_after;
            PELA(j, :) = logic*1;
        end
        sumPELA = sum(PELA, 1);
        proportion(i,:) = sumPELA/size(earthquake_bubble, 1);
        event_number(i) = earthquake_number;
    end
    
    p = proportion;
    q = 1-p;
    sqrtn = sqrt(event_number);
    denominator = p.*q;
    numerator = sqrtn.*abs((proportion-0.5));
    
    Z = numerator./denominator;
    
    index = nan(16, 45);
    index2 = nan(16, 45);
    for i = 1: 16
        
        extP = p(i, :);
        extZ = Z(i, :);
        
        u1 = find(extZ > 2.58 & extP > 0.5);
        u2 = find(isinf(extZ) & extP > 0.5);
        u = union(u1, u2);
        
        u3 = find(extZ > 1.96 & extP > 0.5);
        
        index(i, u) = 1;
        index2(i, u3) = 1;
        
    end
    
    u4 = find(~isnan(index));
    index2(u4) = nan;
    
    fig = figure('position', sz, 'visible', 'off');
    left_color = [.5 .5 0];
    right_color = [1 .5 .5];
    set(fig, 'defaultAxesColorOrder', [left_color; right_color]);
    
    imagesc(1:1:45, (6.5:.1:8)', proportion)
    set(gca, 'ydir', 'normal')
    colormap jet
    hold on
     for i = 1: 16
         plot([0 46], [6.55+0.1*(i-1) 6.55+0.1*(i-1)],...
             'color', [0.5 0.5 0.5], 'linewidth', .1)
     end
     
     for i = 1: 45
         plot([i+0.5 i+0.5], [6 8.5],...
             'color', [0.7 0.2 0.2], 'linewidth', .1)
     end
    
    set(gca, 'ytick', 6.5:0.1:8)
    set(gca, 'fontsize', 20)
    set(gca, 'xtick', 0:3:45)
    
    h = gca;
    xtickangle(0);
    xlabel('Days to Earthquake')
    colorbar
    colormap jet
    ylabel('Magnitude')
    
    yyaxis right
    B = linspace(0,1,33);
    yytk = B(2:2:end);
    set(gca, 'ytick', yytk, 'yticklabel', event_number)
    ylabel('Earthquake Number')
    
    cb = colorbar;
    ylabel(cb, 'Proportion')
    M = 6.5:0.1:8;
    M = M';
    %
    yyaxis left
    hold on
%     P1 = plot(1:1:45, M.*index, 'ko', 'markersize', 6, 'markerface', 'k');
%     P2 = plot(1:1:45, M.*index2, 'ko', 'markersize', 6);
%     legend([P1(1), P2(1)], {'P < 0.01', 'P < 0.05'},...
%         'Position',[0.84 0.027 0.067 0.069])
    
    tenton = 10*(power - floor(power));
    
    title(['ROCSAT-1 \itAi > 10^' num2str(floor(power)) '^.^' num2str(tenton)])
    
    set(gcf,'paperpositionmode','auto')
    print('-dpng','-r300', [output 'Noztest_' num2str(power*10) '_area5.jpg'])
    close all
end

