%============================
% d file transform to o file
%============================
clear;clc

Year = 2011;
path_crx2rnx = 'E:\gnss_tsai\crx2rnx.exe';
% ------------------------------------------------------------------------
y = num2str(Year);
ye2 = y(end-1:end);

for d = 221
    
    cd (['E:\' num2str(Year) '\' num2str(d, '%03d')])
    fn = dir(['*' ye2 'd']);

    for i = 1: size(fn,1)
        filename = fn(i).name;
        st = filename(1:4);
        if exist(filename,'file')
            eval(['!',path_crx2rnx,' ',fn(i).name,' -> ',st,num2str(d,'%03d'),'0.',ye2,'o'])
        end
%         pause(10)
        delete(filename)
        disp([num2str(d, '%03d') '_' filename])
    end   
end