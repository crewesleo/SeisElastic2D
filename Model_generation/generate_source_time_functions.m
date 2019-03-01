clear;
clc

dt=5e-4;
fdom=12;
tlength=0.9995;

[w,tw]=wavemin(dt,fdom,tlength);
%w(1001:end,1)=0;
w=w./max(abs(w));
figure(1);
plot(tw,w);hold on;


fid=fopen('source_wavelet.dat','w');
for ii=1:length(w)
fprintf(fid,'%3.6f	%3.6f',tw(ii,1), w(ii,1));
fprintf(fid,'\n');
end
fclose(fid);






