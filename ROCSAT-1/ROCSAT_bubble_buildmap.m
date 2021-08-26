clear
clc

addpath('E:\Matlabtoolbox')
input = 'E:\Satellite\ROCSAT-1\analyzed\';
output = 'E:\OneDrive - 國立中央大學\myPaper\paper1\JGR\ROCSAT-1\MAP\';
coordinate = 'Magnetic';

TimePeriod = [20, 24];

LON = -179.5 : 5 : 179.5;
LAT = 64.5 : -2.5 : -64.5;
LAT2 = 90 : -2.5 : -90;
lonr = 10;
latr = 2;

Month = {char('January'), char('February'), char('March'), char('April'),...
    char('May'), char('June'), char('July'), char('August'),...
    char('September'), char('October'), char('November'), char('December')};

for k = 1 : 12
   
    fn = [];
    for year = 1999 : 2004
        switch mod(year, 4)
            case {1, 2, 3}
                January = 1 : 31;
                February = 32 : 59;
                March = 60 : 90;
                April = 91 : 120;
                May = 121 : 151;
                June = 152 : 181;
                July = 182 : 212;
                August = 213 : 243;
                September = 244 : 273;
                October = 274 : 304;
                November = 305 : 334;
                December = 335 : 365;
            case 0
                January = 1 : 31;
                February = 32 : 60;
                March = 61 : 91;
                April = 92 : 121;
                May = 122 : 152;
                June = 153 : 182;
                July = 183 : 213;
                August = 214 : 244;
                September = 245 : 274;
                October = 275 : 305;
                November = 306 : 335;
                December = 336 : 366;
        end
        s = Month{k};
        a = numel(eval(s));
        b = nan(a, 1);
        b(:, :) = year;
        b = [b, eval(s)'];
        fn = [fn; num2str(b(:, 1)), num2str(b(:, 2), '%03d')];
    end

    GLON = cell(size(fn,1), 1);
    GLAT = cell(size(fn,1), 1);
    HHTpower = cell(size(fn,1), 1);
    Time = cell(size(fn,1), 1);
    Noi = cell(size(fn,1), 1);
    STD = cell(size(fn,1), 1);
    deltaN = cell(size(fn,1), 1);
    year = cell(size(fn,1), 1);
    month = cell(size(fn,1), 1);
    day = cell(size(fn,1), 1);
    
    for i = 1 : length(fn)
        
        day_data = [];
        f = dir([input fn(i, :) '*mat']);
        if isempty(f);continue;end
        load([input f.name])
        filename = f.name;
        disp(f.name)
        UT = data.Time/3600;
        lon = data.GLON;
        LT = UT + data.GLON/15;
        LT(LT>24) = LT(LT>24) - 24;
        LT(LT<0) = LT(LT<0) + 24;

        u = find(LT >= TimePeriod(1) & LT < TimePeriod(2));
        lon(lon>180) = lon(lon>180) - 360;
        
        GLON{i} = lon(u);
        GLAT{i} = data.GLAT(u);
        HHTpower{i} = data.htpower(u);
        Time{i} = data.Time(u);
        deltaN{i} = linearscale.deltaNi(u);
        Noi{i} = linearscale.Noi(u);
        STD{i} = linearscale.STD(u);
        
        DataLength = length(lon(u));
        y = str2double(filename(1 : 4));
        doy = str2double(filename(5 : 7));
        date = datetime(y, 1, doy);
        m = date.Month;
        d = date.Day;
                
        year{i} = ones(DataLength, 1)*y;
        month{i} = ones(DataLength, 1)*m;
        day{i} = ones(DataLength, 1)*d;
    end
    
    GLON = cell2mat(GLON);
    GLAT = cell2mat(GLAT);
    HHTpower = cell2mat(HHTpower);
    Time = cell2mat(Time);
    deltaN = cell2mat(deltaN);
    Noi = cell2mat(Noi);
    STD = cell2mat(STD);
    year = cell2mat(year);
    month = cell2mat(month);
    day = cell2mat(day);
    
    switch coordinate
        case 'Magnetic'
            cd('E:\Matlabtoolbox\Geo2mag');
            H = GEO2MAG_APX(year, month, day, GLAT, GLON);
            u = find(H(:, 4) < 15 & H(:, 4) > -15);
            latdata = H(u, 4);
        case 'Geographic'
            u = find(GLAT < 35 & GLAT > -35);
            latdata = GLAT(u);
    end

    londata = GLON(u);
    Ai = sqrt(HHTpower(u));
    ambient = Noi(u);
    stddev = STD(u);
    deltaN = deltaN(u);
    
    Ai_MAP = mymap(LON, LAT, lonr, latr, londata, latdata, Ai, 'median');
    Noi_MAP = mymap(LON, LAT, lonr, latr, londata, latdata, ambient, 'median');
    STD_MAP = mymap(LON, LAT, lonr, latr, londata, latdata, stddev, 'median');
    Abs_MAP = mymap(LON, LAT, lonr, latr, londata, latdata, deltaN, 'median');

    if TimePeriod(1) == 0
        name = ['ROC_bbMap_PostMidnight_' num2str(k, '%02d')];
    else
        name = ['ROC_bbMap_PreMidnight_' num2str(k, '%02d')];
    end
    
    save([output name], 'Ai_MAP', 'Noi_MAP', 'STD_MAP', 'Abs_MAP')
end
