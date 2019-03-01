clear;
clc;

fid=fopen('sources.dat','w');

ns=8;
s_dx=130;
xs0=170;
xs=xs0:s_dx:(ns-1)*s_dx+xs0;
s_index=100:ns-1+100;

zs=610;

line1=[0 0 0];
for ii=1:ns
fprintf(fid,'%3.6f	%3.6f',zs,xs(ii));
fprintf(fid,'\n');
end

ns=8;
s_dx=130;
xs0=170;
xs=xs0:s_dx:(ns-1)*s_dx+xs0;
s_index=100:ns-1+100;

zs=15;

line1=[0 0 0];
for ii=1:ns
fprintf(fid,'%3.6f	%3.6f',zs,xs(ii));
fprintf(fid,'\n');
end
fclose(fid);

