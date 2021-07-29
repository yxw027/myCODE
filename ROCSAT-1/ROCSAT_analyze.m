clear
clc

addpath('E:\Matlabtoolbox')
input = 'E:\Satellite\ROCSAT-1\mat\';
output = 'E:\Satellite\ROCSAT-1\analyzed\';
warning('off')
fn = dir([input '*mat']);

for i = 1: 1200
    load([input fn(i).name])
    filename = fn(i).name;
    disp(filename)
    
    Ni = data.LogN;
    Ni = 10.^Ni;
    
    [ed, hd, ion_dens]= delnan(Ni);
    sm_ion = smooth_orbitdata(ion_dens, 1);
    
    [sigma_log, Noi_log, STD_log] = sigmaindex(sm_ion, 'log');
    [sigma_linear, Noi_linear, STD_linear] = sigmaindex(sm_ion, 'linear');
    
    htpower = myHHT(sm_ion);
    sm_ion = [nan(hd,1); sm_ion; nan(ed,1)];
    sigma_log = [nan(hd,1); sigma_log; nan(ed,1)];
    Noi_log = [nan(hd,1); Noi_log; nan(ed,1)];
    STD_log = [nan(hd,1); STD_log; nan(ed,1)];
    
    sigma_linear = [nan(hd,1); sigma_linear; nan(ed,1)];
    Noi_linear = [nan(hd,1); Noi_linear; nan(ed,1)];
    STD_linear = [nan(hd,1); STD_linear; nan(ed,1)];
%--------------------------------------------------------------------------
    Nbar = smooth(sm_ion, 10, 'moving');
    dNi = sm_ion - Nbar;
    count = 0;
    deltaNi = nan(size(sm_ion,1),1);
    bkgd = nan(size(sm_ion,1),1); % background density
    for k = 5: length(sm_ion)-4
        count = count + 1;
        ni = dNi(k-4: k+4);
        if ~isempty(find(isnan(ni),1))
            deltaNi(k,1) = nan;
            bkgd(k,1)  = nan;
            continue
        end
        deltaNi(k,1) = sqrt(nanmean(ni.^2));
        bkgd(k,1)  = nanmean(Nbar(k-4: k+4));
    end
%--------------------------------------------------------------------------
    Nbar = smooth(sm_ion, 10, 'moving');
    dNi_log = log10(sm_ion) - log10(Nbar);
    count = 0;
    deltaNi_log = nan(size(sm_ion,1),1);
    bkgd = nan(size(sm_ion,1),1); % background density
    for k = 5: length(sm_ion)-4
        count = count + 1;
        ni = dNi_log(k-4: k+4);
        if ~isempty(find(isnan(ni),1))
            deltaNi_log(k,1) = nan;
            bkgd(k,1)  = nan;
            continue
        end
        deltaNi_log(k,1) = sqrt(nanmean(ni.^2));
    end
%--------------------------------------------------------------------------
    
    deltaNi = [nan(hd,1); deltaNi; nan(ed,1)];
    deltaNi_log = [nan(hd,1); deltaNi_log; nan(ed,1)];
    bkgd = [nan(hd,1); bkgd; nan(ed,1)];
    
    data.htpower = htpower;
    data.sm_ion = sm_ion;
    data.Time = (1:1:86400)';
    data.bkgd = bkgd;
    
    linearscale.sigma = sigma_linear;
    linearscale.Noi = Noi_linear;
    linearscale.STD = STD_linear;
    linearscale.deltaNi = deltaNi;
    
    logscale.sigma = sigma_log;
    logscale.Noi = Noi_log;
    logscale.STD = STD_log;
    logscale.deltaNi = deltaNi_log;
    
    save([output filename '.mat'], 'data', 'linearscale', 'logscale')
    disp([filename ' profile is done.'])
end