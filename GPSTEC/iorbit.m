function iorbit(Y, D, rate, OUTFILE1)

ORBPATH = 'ORB/';

if ~exist([ORBPATH OUTFILE1],'file')
    INFILE2 = OUTFILE1([1:8 12:end]); % CODwwwwd.mat
    if exist([ORBPATH INFILE2],'file')
        sp3cubic([ORBPATH INFILE2],rate); % Create CODwwwwdiRR.mat
    elseif exist([ORBPATH INFILE2(1:8) '.EPH'],'file') % CODwwwwd.EPH
        sp32mat(w,wd); % data CODwwwwd.mat
        sp3cubic([ORBPATH INFILE2],rate); % Create CODwwwwdiRR.mat
    else
        for i=[115 114 117] % ASCII numbers of 's', 'r' and 'u'
            OUTFILE1 = ['ig' char(i) W WD 'i' R '.mat']; % ig[sru]wwwwdiRR.mat
            if exist([ORBPATH OUTFILE1],'file')
                break;
            else
                INFILE2 = OUTFILE1([1:8 12:end]); % ig[sru]wwwwd.mat
                if exist([ORBPATH INFILE2],'file')
                    sp3cubic([ORBPATH INFILE2],rate); % Create ig[sru]wwwwdiRR.mat
                    break;
                elseif exist([ORBPATH INFILE2(1:8) '.sp3'],'file') % ig[sru]wwwwd.sp3
                    sp32mat(w,wd); % Create ig[sru]wwwwd.mat
                    sp3cubic([ORBPATH INFILE2],rate); % Create ig[sru]wwwwdiRR.mat
                    break;
                end
            end
        end
    end
end