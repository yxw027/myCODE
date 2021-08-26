function getgnssdata(yr,doy, MISSION)

MISSION = upper(MISSION);
codeftp  = 'ftp://ftp.aiub.unibe.ch/CODE/';
brdcftp = 'ftp://cddis.gsfc.nasa.gov/glonass/data/daily/';
galftp = 'ftp://ftp.gfz-potsdam.de/GNSS/products/mgex/';
[yr,mon,day] = doy2ymd(yr,doy);
Y = num2str(yr,'%04u'); Y2 = Y(3:4);
D = num2str(doy,'%03u'); % Day of year in 3-char
M = num2str(mon,'%02u');
% DAY = num2str(day,'%02u');
[w,wd] = ymd2gps(yr,mon,day);
W = num2str(w,'%04u');
WD = int2str(wd);
rate = 30;

ORBPATH = 'ORB/';

if ~exist(ORBPATH, 'dir')
   mkdir(ORBPATH);
end

EPMPATH = 'EPM/';
if ~exist(EPMPATH, 'dir')
   mkdir(EPMPATH);
end

DCBPATH = 'DCB/';
if ~exist(DCBPATH, 'dir')
   mkdir(DCBPATH);
end
% orbfile
switch MISSION
    case {'GPS','GLONASS'}
        orbfile = ['COD' W WD '.EPH.Z'];
        if ~exist([ORBPATH orbfile],'file') && ~exist([ORBPATH orbfile(1:end-2)],'file')
            system(['wget -q ' codeftp  Y '/' orbfile]);
            system(['gzip -d -f ' orbfile ]);
            movefile(orbfile(1:end-2), ORBPATH);
            sp32mat(w,wd);
            sp3cubic([ORBPATH 'COD' W WD '.mat'], rate);
        end
     case {'GALILEO'}
         orbfile = ['gbm' W WD '.sp3.Z'];
        if ~exist([ORBPATH orbfile],'file') && ~exist([ORBPATH orbfile(1:end-2)],'file')
            system(['wget -q ' galftp  W '/' orbfile]);
            system(['gzip -d -f ' orbfile ]);
            movefile(orbfile(1:end-2),ORBPATH);
            sp32mat_e(w,wd);
            sp3cubic([ORBPATH 'gbm' W WD '.mat'], rate);
         end
end

% DCB file
switch MISSION
       case {'GPS','GLONASS'}
            dcbfile = ['P1P2' Y2 M '_ALL.DCB']; % P1P2YYMM_ALL.mat
            INFILE6 = [DCBPATH dcbfile]; 
            if ~exist(INFILE6,'file')
                dcbfile1 = [dcbfile '.Z'];
                system(['wget -q ' codeftp  Y '/' dcbfile1]);
                system(['gzip -d -f ' dcbfile1 ]);
                movefile(dcbfile, DCBPATH);
                dcb2mat([DCBPATH dcbfile]);
            end
            dcbfile1 = ['P1C1' Y2 M '_RINEX.DCB']; % P1C1YYMM_RINEX.mat
            INFILE7 = [DCBPATH dcbfile1];
            if ~exist(INFILE7,'file')
              %  warning([INFILE7 ' not found. Downloading now!!']);
                dcbfile = [dcbfile1 '.Z'];
                system(['wget -q ' codeftp  Y '/' dcbfile]);
                system(['gzip -d -f ' dcbfile ]);
                movefile(dcbfile1, DCBPATH);
                dcb2mat([DCBPATH dcbfile1]);
            end
end

% broadcast ephemeris files 
if strcmp(MISSION,'GLONASS')
    brdcfile = ['brdc' D '0.' Y2 'g'];
    if ~exist(brdcfile, 'file')
        system(['wget -q ' brdcftp  Y '/' D '/' Y2 'g/' brdcfile '.Z']);
        system(['gzip -d -f ' brdcfile '.Z' ]);
        movefile(brdcfile,EPMPATH);
    end
end




