function [clight, f1, f2, ionfac, re, h1, h2, nsat, minobs, km2m] = readparm(str)

fileID = fopen(str);
C = textscan(fileID,'%s %f');
fclose(fileID);
for i = 1: 10
eval([C{1}{i} '=' num2str(C{2}(i)) ';'])
end