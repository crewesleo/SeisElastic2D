clear;
clc;

ctypes={'div', 'seq', 'qual'};
cnames{1,:}={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};
cnames{2,:}={'Blues','BuGn','BuPu','GnBu','Greens','Greys','Oranges','OrRd','PuBu','PuBuGn','PuRd',...
             'Purples','RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd'};
cnames{3,:}={'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'};cnames{1,:}={'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'};


load ./ModelParams/Vp_true_small.mat;
load ./ModelParams/Vs_true_small.mat;
load ./ModelParams/Density_true_small.mat;
load ./ModelParams/Epsilon_true_small.mat;
load ./ModelParams/Delta_true_small.mat;
load ./ModelParams/Theta_true_small.mat;

Vp=Vp_true_small(1:2:end,750:3:1125);
Vs=Vs_true_small(1:2:end,750:3:1125);
Rho=Density_true_small(1:2:end,750:3:1125);
epsilon=Epsilon_true_small(1:2:end,750:3:1125);
delta=Delta_true_small(1:2:end,750:3:1125);
theta=Theta_true_small(1:2:end,750:3:1125);
theta=theta/2;
[nz,nx]=size(Vp);

Rho=gaussian_smoother(Rho,1:nx,1:nz,10);

% unroatated elastic constants
C11U=Rho.*Vp.^2.*(1+2.*epsilon);
C13U=sqrt(2*Rho.*(Vp.^2).*delta.*(Rho.*Vp.^2-Rho.*Vs.^2)+(Rho.*Vp.^2-Rho.*Vs.^2).^2)-Rho.*Vs.^2;
C33U=Rho.*Vp.^2;
C55U=Rho.*Vs.^2;
THETA=theta*pi/180;

p_1=C11U-2*C13U+C33U-4*C55U;
p_2=C11U-C33U;

%
C11=(3*C11U+2*C13U+3*C33U+4*C55U+4*p_2.*cos(2*THETA)+p_1.*cos(4*THETA))/8;
C13=(C11U+6*C13U+C33U-4*C55U-p_1.*cos(4*THETA))/8;
C15=(p_2.*sin(2*THETA)+p_1.*sin(4*THETA))/4;
C33=(3*C11U+2*C13U+3*C33U+4*C55U-4*p_2.*cos(2*THETA)+p_1.*cos(4*THETA))/8;
C35=(p_2.*sin(2*THETA)-p_1.*sin(4*THETA))/4;
C55=(C11U-2*C13U+C33U+4*C55U-p_1.*cos(4*THETA))/8;
C12=ones(nz,nx)*0;
C23=ones(nz,nx)*0;
C25=ones(nz,nx)*0;

xx=1:nx;
zz=1:nz;
xx=xx.*5;
zz=zz.*5;

figure;
subplot(1,6,1)
imagesc(xx,zz,Vp)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,2)
imagesc(xx,zz,Vs)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,3)
imagesc(xx,zz,Rho)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,4)
imagesc(xx,zz,epsilon)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,5)
imagesc(xx,zz,delta)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,6)
imagesc(xx,zz,theta)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
colors = colormap(cbrewer(ctypes{1},cnames{1}{5}, 100));
% colormap('jet')

figure;
subplot(1,6,1)
imagesc(xx,zz,C11)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,2)
imagesc(xx,zz,C13)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,3)
imagesc(xx,zz,C15)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,4)
imagesc(xx,zz,C33)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,5)
imagesc(xx,zz,C35)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,6,6)
imagesc(xx,zz,C55)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
colors = colormap(cbrewer(ctypes{1},cnames{1}{5}, 100));
% colormap('jet')

figure;
subplot(1,5,1)
imagesc(xx,zz,C11U)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,5,2)
imagesc(xx,zz,C13U)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,5,3)
imagesc(xx,zz,C33U)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,5,4)
imagesc(xx,zz,C55U)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
subplot(1,5,5)
imagesc(xx,zz,THETA)
xlabel('Distance (m)','fontname','arial');
ylabel('Depth (m)','fontname','arial');
colorbar
colors = colormap(cbrewer(ctypes{1},cnames{1}{5}, 100));


%%
% colors = colormap(cbrewer(ctypes{1},cnames{1}{4}, 100));

[nz,nx]=size(Vp);

dx=5;
dz=5;

x0=0;
z0=0;

xend=(nx-1)*dx;
zend=(nz-1)*dz;

min_vp=min(min(Vp));
max_vp=max(max(Vp));

min_vs=min(min(Vs));
max_vs=max(max(Vs));

min_rho=min(min(Rho));
max_rho=max(max(Rho));


min_c11=min(min(C11));
max_c11=max(max(C11));

min_c13=min(min(C13));
max_c13=max(max(C13));

min_c15=min(min(C15));
max_c15=max(max(C15));

min_c33=min(min(C33));
max_c33=max(max(C33));

min_c35=min(min(C35));
max_c35=max(max(C35));

min_c55=min(min(C55));
max_c55=max(max(C55));

min_c12=min(min(C12));
max_c12=max(max(C12));

min_c23=min(min(C23));
max_c23=max(max(C23));

min_c25=min(min(C25));
max_c25=max(max(C25));


Vp=Vp';
Vs=Vs';
Rho=Rho';

C11=C11';
C13=C13';
C15=C15';
C33=C33';
C35=C35';
C55=C55';
C12=C12';
C23=C23';
C25=C25';

vp=Vp(:);
vs=Vs(:);
rho=Rho(:);

c11=C11(:);
c13=C13(:);
c15=C15(:);

c33=C33(:);
c35=C35(:);
c55=C55(:);

c12=C12(:);
c23=C23(:);
c25=C25(:);

column6=ones(size(vp))*9998;
column7=ones(size(vp))*9998;


xx=0:dx:xend;
x_coordinate=ones(size(Vp'));
for iii=1:nz
    x_coordinate(iii,:)=xx';
end
x_coordinate=x_coordinate';
x_coordinate=x_coordinate(:);


for iii=1:nz
    z_coordinate(iii,:)=ones(size(xx))*(iii-1)*dz;
end

z_coordinate=z_coordinate';
z_coordinate=z_coordinate(:);


model_file=zeros(length(x_coordinate),16);

model_file(:,1)=x_coordinate;
model_file(:,2)=z_coordinate;
model_file(:,3)=vp;
model_file(:,4)=vs;
model_file(:,5)=rho;
model_file(:,6)=column6;
model_file(:,7)=column7;

model_file(:,8)=c11;
model_file(:,9)=c13;
model_file(:,10)=c15;

model_file(:,11)=c33;
model_file(:,12)=c35;
model_file(:,13)=c55;

model_file(:,14)=c12;
model_file(:,15)=c23;
model_file(:,16)=c25;

% save model_file_r.xyz -ascii model_file

fid=fopen('model_file_bp_tti_true.xyz','w');
[m,n]=size(model_file);
line1=[x0 z0 xend zend];
line2=[dx dz];
line3=[nx nz];
line4=[min_vp max_vp min_vs max_vs min_rho max_rho 9998 9998 9998 9998 min_c11 max_c11 min_c13 max_c13 min_c15 max_c15 min_c33 max_c33 min_c35 max_c35 min_c55 max_c55 min_c12 max_c12 min_c23 max_c23 min_c25 max_c25];

for ii=1:length(line1)
fprintf(fid,'%3.1f ',line1(ii));
end
fprintf(fid,'\n');

for ii=1:length(line2)
fprintf(fid,'%3.1f ',line2(ii));
end
fprintf(fid,'\n');

for ii=1:length(line3)
fprintf(fid,'%3.1f ',line3(ii));
end
fprintf(fid,'\n');

for ii=1:length(line4)
fprintf(fid,'%3.1f ',line4(ii));
end
fprintf(fid,'\n');

for i=1:m
    for j=1:n

        fprintf(fid,'%3.1f ',model_file(i,j));
        if j==n
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);
