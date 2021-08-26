clear;clc

% cd('E:\gnss_tsai')
% addpath('E:\gnss_tsai');

Year = 2011;
% -------------------------------------------------------------------------
for doy = 221
    
    input = ['E:\' num2str(Year) '\' num2str(doy, '%03d') '\']; % o-file directory
    output = ['E:\TEC30\' num2str(Year) '\' num2str(doy, '%03d') '\'];
    if ~exist(output, 'dir')
        mkdir(output);
    end
    
    hrng = [0 24]; % hour range
    rate = 30; % sampling rate
    minel = 30; %  minimum (cutoff) elevation angle of satellites
    ippalt = 450; % ionospheric altitude
    Y = num2str(Year,'%04d');
    mission ={'GPS', 'GLONASS', 'Galileo'};
    m1 = 1;% mission 1
    m2 = 1;% mission 2
    for mm = m1: m2
        MISSION = mission{mm};
        switch MISSION
            case 'GPS'
                fun = @gpstec;
            case 'GLONASS'
                fun = @glotec;
            case 'Galileo'
                fun = @galtec;
        end
        
        D = num2str(doy,'%03d');
        getgnssdata(Year, doy, MISSION); % download necessary data
        fn = dir([input '*0.' Y(3:4) 'o']);
        for i = 1 : size(fn,1)
            infile = fn(i).name;
            Site = infile(1:4);
            if exist([output infile(1:7) '.mat'], 'file')
                continue;
            end
            try
                fun(Site, Year, doy, hrng, rate, minel, ippalt, input, output);
            catch
                disp([Site ' wrong !!']);
                continue;
            end
        end
    end
end