function [OrbitTime, OrbitLat, OrbitLon, OrbitAlt, OrbitData1,...
    OrbitData2, OrbitData3] = ExtractandDistinguish(Time, Lat, Lon,...
    Alt, data1, data2, data3, epicenterLon, epicenterLat, UT, minUT, maxUT)

% ExtractandDistinguish(x) returns extracted and distingiushed data of DMSP
%   AUTHOR	:	Chun-Yen Huang
%   SINCE	:	2021/8/6
%   VERSION	:	1.0 2021/8/6

u1 = find(UT > minUT & UT < maxUT);
u2 = find(Lon > epicenterLon - 20 & Lon < epicenterLon + 50);
u = intersect(u1, u2);

Time = Time(u, :);
Lat = Lat(u, :);
Lon = Lon(u, :);
Alt = Alt(u, :);
data1 = data1(u, :);
data2 = data2(u, :);
data3 = data3(u, :);

U = [];
for i = 1: length(Lat)-1
    a = abs(Lat(i+1)-Lat(i));
    if a > 10
        U = [U,i];
    else
    end
end

OrbitTime = cell(length(U)+1, 1);
OrbitLat = cell(length(U)+1, 1);
OrbitLon = cell(length(U)+1, 1);
OrbitAlt = cell(length(U)+1, 1);
OrbitData1 = cell(length(U)+1, 1);
OrbitData2 = cell(length(U)+1, 1);
OrbitData3 = cell(length(U)+1, 1);
for i = 1: size(OrbitLat,1)
    if i == 1 && ~isempty(U)
        OrbitTime{i, 1} = Time(1:U(i),:);
        OrbitLat{i, 1} = Lat(1:U(i),:);
        OrbitLon{i, 1} = Lon(1:U(i),:);
        OrbitAlt{i, 1} = Alt(1:U(i),:);
        OrbitData1{i, 1} = data1(1:U(i),:);
        OrbitData2{i, 1} = data2(1:U(i),:);
        OrbitData3{i, 1} = data3(1:U(i),:);
    elseif i > 1 && i < size(OrbitData1,1)
        OrbitTime{i, 1} = Time(U(i-1)+1:U(i),:);
        OrbitLat{i, 1} = Lat(U(i-1)+1:U(i),:);
        OrbitLon{i, 1} = Lon(U(i-1)+1:U(i),:);
        OrbitAlt{i, 1} = Alt(U(i-1)+1:U(i),:);
        OrbitData1{i, 1} = data1(U(i-1)+1:U(i),:);
        OrbitData2{i, 1} = data2(U(i-1)+1:U(i),:);
        OrbitData3{i, 1} = data3(U(i-1)+1:U(i),:);
    elseif isempty(U)
        OrbitTime{i, 1} = Time;
        OrbitLat{i, 1} = Lat;
        OrbitLon{i, 1} = Lon;
        OrbitAlt{i, 1} = Alt;
        OrbitData1{i, 1} = data1;
        OrbitData2{i, 1} = data2;
        OrbitData3{i, 1} = data3;
    else
        OrbitTime{i, 1} = Time(U(i-1)+1:end,:);
        OrbitLat{i, 1} = Lat(U(i-1)+1:end,:);
        OrbitLon{i, 1} = Lon(U(i-1)+1:end,:);
        OrbitAlt{i, 1} = Alt(U(i-1)+1:end,:);
        OrbitData1{i, 1} = data1(U(i-1)+1:end,:);
        OrbitData2{i, 1} = data2(U(i-1)+1:end,:);
        OrbitData3{i, 1} = data3(U(i-1)+1:end,:);
    end
end