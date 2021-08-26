function [dflag] = DLOBS(SD, Y, D, rate, inpath)
Y2 = Y(3:4);
R = num2str(rate,'%02u');
obsmatdir =['OBS/gps/' Y '.' D '/'];
% Deal with the O-file(s)
OUTFILE = [SD '0.mat']; % SSSSDDD0.mat
switch rate
    case {30}
        INFILE = [inpath OUTFILE(1:7) '0.' Y(3:4) 'o']; % SSSSDDD0.YYo
    case {15}
        INFILE = [inpath OUTFILE(1:7) '1.' Y(3:4) 'o']; % SSSSDDD1.YYo
    case {1}
        INFILE = [inpath OUTFILE(1:7) '0.' Y(3:4) 'o']; % SSSSDDD2.YYo
    case {0.01}
        INFILE = [inpath OUTFILE(1:7) '2.' Y(3:4) 'o']; % SSSSDDD2.YYo
end

if exist(INFILE,'file')
    rv = floor(str2double(textread(INFILE, '%s',1)));
    switch rv
        case {2}
            tmpfile = [OUTFILE(1:9) Y2 'o.t1'];
            system(['teqc -phc -O.obs L1L2C1P1P2 -O.dec ' R 's ' INFILE ' > ' tmpfile]); % SSSSDDD0.YYo to SSSSDDD0.YYo.t1
            sflag = obs2mat(rv, tmpfile); % SSSSDDD0.YYo.t1 to SSSSDDD0.mat
            if sflag
                delete(tmpfile);
                warning('The site is in the skip list or no any observation data found.');
                return;
            end
            delete(tmpfile); % Remove the temporary file
            movefile(OUTFILE,obsmatdir);
        case 3
            sflag = obs2mat(INFILE); % SSSSDDD0.YYo to SSSSDDD0.mat
            if sflag
                warning('The site is in the skip list or no observation data found.');
                return;
            end
            movefile(OUTFILE,obsmatdir);
        otherwise
            error('The RINEX version of the input file should be 2.00 or later.');
    end
%     if rate == 30
        dflag = true;
%     else
%         dflag = false;
%     end
else
    error(['Cannot find ' INFILE '.']);
end
