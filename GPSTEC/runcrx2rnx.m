


yr =2016;
doy=249;
Y = num2str(yr,'%04d');
D = num2str(doy,'%03d');
inpath =['/nishome/michelle/RINEX/' Y '.' D '/'];
f= dir([inpath '*.Z']);
for i=1:length(f)
    infile = f(i).name;
    Site = infile(1:4);
	crx2rnx(Site,yr,doy,inpath);

end
