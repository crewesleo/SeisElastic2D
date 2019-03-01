clear;
clc;

fid=fopen('STATIONS','w');

ns=81;
s_dx=15;
xs0=25;
xs=xs0:s_dx:(ns-1)*s_dx+xs0;
s_index=100:ns-1+100;

line1=[0 0];
for ii=1:ns
fprintf(fid,'%s %s  %3.2f %3.2f ',strcat('S0',num2str(s_index(ii))),'AA',15,xs(ii));
fprintf(fid,'%3.2f ',line1);

fprintf(fid,'\n');
end

ns=81;
s_dx=-15;
xs0=1225;
xs=xs0:s_dx:(ns-1)*s_dx+xs0;
s_index=100:ns-1+100;

line1=[0 0];
for ii=1:ns
fprintf(fid,'%s %s  %3.2f %3.2f ',strcat('S0',num2str(s_index(ii))),'AA',610,xs(ii));
fprintf(fid,'%3.2f ',line1);

fprintf(fid,'\n');
end

fclose(fid);