clear
clc

input = 'E:\Satellite\ROCSAT-1\raw\';
output = 'E:\Satellite\ROCSAT-1\mat\';

for year = 1999: 2004
    cd ([input num2str(year)])
    fn = dir('*.dat');
    
    for i = 1: size(fn,1)
        disp([num2str(year) '_' num2str(i)])
        temp = importdata([input num2str(year) '\' fn(i).name],' ',1);
        filename = fn(i).name;
        doy = str2double(filename(3:5));
        time = temp.textdata(2:end,2);
        time = split(time, ':');
        A = cellfun(@str2double, time, 'UniformOutput', false);
        A = cell2mat(A);
        
        L = size(time,1);
        Date = datetime([year*ones(L,1), ones(L,1), doy*ones(L,1), A]);
        Date = datevec(Date);
        second_of_day = A(:,1)*3600 + A(:,2)*60 + A(:,3) + 1;
        
        header = regexp(temp.textdata(1), '\s+', 'split');
        header = header{1};
        header{16} = 'O';
        tempdata = temp.data;
        
        for j = 1: 16
            var = header{j+7};
            eval([var '=nan(86400, 1);'])
            datas = tempdata(:, j);
            eval([var '(second_of_day(second_of_day>=1 & second_of_day<=86400)) = datas(second_of_day>=1 & second_of_day<=86400);'])
            eval(['data.' var ' = ' var ';'])
        end
        
        
        Time = nan(86400, 6);
        Time(second_of_day(second_of_day>=1 & second_of_day<=86400), :) = Date(second_of_day>=1 & second_of_day<=86400, :);
        
        save([output num2str(year) filename(3:5)], 'data', 'Time', 'header')
    end
end