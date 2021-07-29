%%
clear
clc

addpath('E:\AOGS_study')
a = load('eqlist_DEMETER_20210711_Global3');
earthquake = a.data2;

input = 'E:\AOGS_study\DEMETER\global_event\';
output = 'E:\AOGS_study\DEMETER\midlow\PEIA_area\';
load([input 'DEMETER_data.mat'])
u = find(GLAT > -35 & GLAT < 35);

GLAT = GLAT(u);
GLON = GLON(u);
HHT = HHT(u);
time = time(u,:);
Ni = Ni(u);
sigma = sigma(u);


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
    eval(['DEMETER_epdata = cell(size(' A ',1), 90);'])
    for j = 1: lght(i)
        disp(j)
        eval(['date = ' A '(j, 1: 3);'])
        
        eval(['eplat = ' A '(j,7);']);
        eval(['eplon = ' A '(j,8);']);
        eval(['epM = ' A '(j,10);']);
        
        doy = day(datetime(date), 'dayofyear');
        eval(['yr = ' A '(j,1);'])
        
        count = 0;
        for k = -45: 1: 45
            count = count + 1;
            disp([num2str(j) '_' num2str(k)])
            
            tempday = datevec(datetime(date)+days(k));
            
            u = find(time(:,1) == tempday(1) & time(:,2) == tempday(2) &...
                time(:,3) == tempday(3));
            
            lat = GLAT(u);
            lon = GLON(u);
            N = Ni(u);
            Si = sigma(u);
            htpower = HHT(u);
            
            R_strain = 10.^(0.43*epM);
            
            uu = find(lat >= eplat - R_strain/110 & lat <= eplat + R_strain/110 &...
                lon >= eplon - R_strain/105 & lon <= eplon + R_strain/105);
            
            DEMETER_epdata{j, count} = [lat(uu), lon(uu), N(uu), htpower(uu), Si(uu)];
        end
    end
    save([output A ],  A , 'DEMETER_epdata')
end

%%
clear
clc

input = 'E:\AOGS_study\DEMETER\global_event\Global_PEIA_range01\';
fn = dir([input '*mat']);

sz = get(0, 'screensize');

datastr = {'data65', 'data66', 'data67', 'data68', 'data69', 'data70', 'data71',...
    'data72', 'data73', 'data74', 'data75', 'data76', 'data77', 'data78'....
    , 'data79', 'data80'};

for i = 14
    
    disp(datastr{i});
    
    load([input fn(i).name])
    
    earthquake_bubble = nan(size(DEMETER_epdata));
    earthquake_number = size(DEMETER_epdata, 1);
    for j = 1: earthquake_number
        for k = 1: 91
            
            a = DEMETER_epdata{j, k};
            if isempty(a); continue; end
            u = find(sqrt(a(:,4))>10^4);
            
            earthquake_bubble(j, k) = numel(u);
            
        end
    end
    earthquake_bubble(isnan(earthquake_bubble)) = 0;
end

[pst, index] = sort(data78(:, 10));
data78 = data78(index, :);


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
    
    lat = num2str(data78(c,7));
    lon = num2str(data78(c,8));
    M = num2str(data78(c, 10));
    strdate = datestr(data78(c,1:6), 26);
    
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
PELA = nan(size(data78,1), 45);
for i = 1: size(data78, 1)
    
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

path_f = 'E:\AOGS_study\DEMETER\';
set(gcf,'paperpositionmode','auto')
print('-dpng','-r300',[path_f 'DEMETER_temporal.jpg'])
close all
    
%%
clear
clc

input = 'E:\AOGS_study\DEMETER\lowlat\PEIA_area\';
output = 'E:\AOGS_study\DEMETER\lowlat\PEIA_area\';

fn = dir([input '*mat']);

    datastr = {'data65', 'data66', 'data67', 'data68', 'data69', 'data70', 'data71',...
        'data72', 'data73', 'data74', 'data75', 'data76', 'data77', 'data78'....
        , 'data79', 'data80'};
sz = get(0, 'screensize');
    
for power = 3.5: 0.1: 5
    proportion = nan(size(fn,1), 45);
    event_number = nan(size(fn,1),1);
    for i = 1: size(fn, 1)
        load([input fn(i).name])
        A = datastr{i};
        
        earthquake_bubble = nan(size(DEMETER_epdata));
        bubble_amplitude = cell(size(DEMETER_epdata));
        earthquake_number = size(DEMETER_epdata, 1);
        for j = 1: earthquake_number
            for k = 1: 91
                
                a = DEMETER_epdata{j, k};
                if isempty(a); continue; end
                u = find(sqrt(a(:,4))/2 >= 10^power);
                
                earthquake_bubble(j, k) = numel(u);
                bubble_amplitude{j, k} = sqrt(a(u,4))/2;
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
%     
    fig = figure('position', sz, 'visible', 'off');
    left_color = [.5 .5 0];
    right_color = [1 .5 .5];
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    
    imagesc(1:1:45, (6.5:.1:8)', proportion)
    caxis([0 1])
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
    P1 = plot(1:1:45, M.*index, 'ko', 'markersize', 6, 'markerface', 'k');
    P2 = plot(1:1:45, M.*index2, 'ko', 'markersize', 6);
    legend([P1(1), P2(1)], {'P < 0.01', 'P < 0.05'},...
        'Position',[0.84 0.027 0.067 0.069])
    
    tenton = 10*(power - floor(power));
    
    title(['DEMETER \itAi > 10^' num2str(floor(power)) '^.^' num2str(tenton)])
    
    saveas(gcf, [output 'ztest_' num2str(power*10)], 'jpg')
    close all
end