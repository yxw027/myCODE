p = gcp;
if isempty(p)
    error('pctexample:backslashbench:poolClosed', ...
        ['This example requires a parallel pool. ' ...
         'Manually start a pool using the parpool command or set ' ...
         'your parallel preferences to automatically start a pool.']);
end
poolSize = p.NumWorkers;

dirpath = ['/nishome/michelle/RINEX/'];
outpath1 = ['/nishome/michelle/TEC/'];
%dirpath = ['/pub/geonet/'];
%dirpath = ['/pub/gnss/Korea_GNSS/'];
yr = 2017;
hrng =[0 24];
rate = 30;
minel = 20;
doy1 = 76;
doy2 = 76;
ippalt = 250;
m1=1;
m2=1;


Y = num2str(yr,'%04d');
mission ={'GPS','GLONASS','Galileo'};
for mm = m1:m2
MISSION = mission{mm};

switch MISSION
    case 'GPS'
       fun = @gpstec;
    case 'GLONASS'
       fun = @glotec;
    case 'Galileo'
       fun = @galtec;
end
%dirpath = ['/pub/gnss/cwb/data/daily/'];
%fid = fopen('wrongsta','w');
%c=1;
for doy = doy1:doy2
    D = num2str(doy,'%03d');
    getgnssdata(yr,doy,MISSION);
    inpath =[dirpath '/' Y '.' D '/' ];
    outpath = [outpath1  MISSION '/' Y '.' D '/'];
    if ~exist(outpath,'dir')
       mkdir(outpath);
    end
    f = dir([inpath '*.*o']);
    lf = length(f);
    parfor (i=1:lf,poolSize)
        infile = f(i).name;
        Site = infile(1:4);
        %if exist([outpath infile(1:7) '.mat'],'file')
        %   continue;
        %end
        try
            %gpstec(Site,yr,doy,hrng,rate,minel);
            fun(Site,yr,doy,hrng,rate,minel,ippalt,inpath,outpath);
        catch
            disp([Site ' wrong !!']);
            %wsta(1,i,c)=1;%fprintf(fid,'station: %s\n',[Site ' ' D ]);
            continue;
        end
    end
    %c = c+1;
end
end
