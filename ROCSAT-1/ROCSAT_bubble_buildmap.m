clear
clc

addpath('E:\Matlabtoolbox')
input = 'E:\Satellite\ROCSAT-1\analyzed\';
output = 'E:\AOGS_study\ROCSAT-1\bubble_MAP\';

LON = -179.5 : 5 : 179.5;
LAT = 64.5 : -2.5 : -64.5;
LAT2 = 90 : -2.5 : -90;
lonr = 10;
latr = 1;

Month = {char('January'), char('February'), char('March'), char('April'),...
    char('May'), char('June'), char('July'), char('August'), char('September'),...
    char('October'), char('November'), char('December')};

for k = 6
   
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
    HHT = cell(size(fn,1), 1);
    sigma = cell(size(fn,1), 1);
    Ni = cell(size(fn,1), 1);
    for i = 1 : length(fn)
        
        day_data = [];
        f = dir([input fn(i, :) '*mat']);
        if isempty(f);continue;end
        load([input f.name])
        disp(f.name)
        UT = data.Time/3600;
        lon = data.GLON;
        LT = UT + data.GLON/15;
        LT(LT>24) = LT(LT>24) - 24;
        LT(LT<0) = LT(LT<0) + 24;

        u = find(LT >= 20 & LT < 24);
        lon(lon>180) = lon(lon>180) - 360;
        
        GLON{i} =lon(u);
        GLAT{i} = data.GLAT(u);
        HHT{i} = data.htpower(u);
        sigma{i} = logscale.sigma(u);
        Ni{i} = data.LogN(u);
    end
    
    GLON = cell2mat(GLON);
    GLAT = cell2mat(GLAT);
    HHT = cell2mat(HHT);
    sigma = cell2mat(sigma);
    Ni = cell2mat(Ni);


    u = find(GLAT < 35 & GLAT > -35);
    latdata = GLAT(u);
    londata = GLON(u);
    Ai = sqrt(HHT(u));
    N = Ni(u);
    Si = sigma(u);
    
    Ai_MAP = mymap(LON, LAT, lonr, latr, londata, latdata, Ai, 'median');
    Ni_MAP = mymap(LON, LAT, lonr, latr, londata, latdata, 10.^N, 'median');
    PAi_MAP = myoccur(LON, LAT, lonr, latr, londata, latdata, Ai, 10^3.7);
    Si_MAP = mysigma(LON, LAT, lonr, latr, londata, latdata, Si, 0.3);

    name = ['F1_HHT_bbMap_' num2str(k, '%02d')];
    
    save([output name], 'Ai_MAP', 'Ni_MAP', 'PAi_MAP', 'Si_MAP')
end
