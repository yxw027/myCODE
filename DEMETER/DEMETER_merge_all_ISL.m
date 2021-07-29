clear
clc

input1 = 'E:\Satellite\DEMETER\ISL_v2\day\';
input2 = 'E:\Satellite\DEMETER\ISL_v2\night\';
output = 'E:\AOGS_study\DEMETER\global_event\';

fn1 = dir([input1 '*mat']);
fn2 = dir([input2 '*mat']);

daylat = cell(size(fn1,1),1);
daylon = cell(size(fn1,1),1);
dayNi = cell(size(fn1,1),1);
dayHHT = cell(size(fn1,1),1);
daysigma = cell(size(fn1,1),1);
daytime = cell(size(fn1,1),1);
for i = 1: size(fn1,1)
    disp(i)
    load([input1 fn1(i).name])
    daylat{i} = tdata.data(:,9);
    daylon{i} = tdata.data(:,10);
    dayNi{i} = rdata.ion;
    dayHHT{i} = linearscale.HHT;
    daysigma{i} = logscale.sigma;
    daytime{i} = tdata.data(:,1:6);
end

daylat = cell2mat(daylat);
daylon = cell2mat(daylon);
dayNi = cell2mat(dayNi);
dayHHT = cell2mat(dayHHT);
daysigma = cell2mat(daysigma);
daytime = cell2mat(daytime);
daysigma = real(daysigma);


nightlat = cell(size(fn2,1),1);
nightlon = cell(size(fn2,1),1);
nightNi = cell(size(fn2,1),1);
nightHHT = cell(size(fn2,1),1);
nightsigma = cell(size(fn2,1),1);
nighttime = cell(size(fn2,1),1);
for i = 1: size(fn2,1)
    disp(i)
    load([input2 fn2(i).name])
    nightlat{i} = tdata.data(:,9);
    nightlon{i} = tdata.data(:,10);
    nightNi{i} = rdata.ion;
    nightHHT{i} = linearscale.HHT;
    nightsigma{i} = logscale.sigma;
    nighttime{i} = tdata.data(:,1:6);
end
nightlat = cell2mat(nightlat);
nightlon = cell2mat(nightlon);
nightNi = cell2mat(nightNi);
nightHHT = cell2mat(nightHHT);
nightsigma = cell2mat(nightsigma);
nighttime = cell2mat(nighttime);

nightsigma = real(nightsigma);

GLAT = [daylat; nightlat];
GLON = [daylon; nightlon];
Ni = [dayNi; nightNi];
HHT = [dayHHT; nightHHT];
sigma = [daysigma; nightsigma];
time = [daytime; nighttime];

savefast([output 'DEMETER_data'], 'time', 'GLAT', 'GLON', 'Ni', 'HHT', 'sigma')
