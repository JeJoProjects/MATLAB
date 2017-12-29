%%                        Assignment - 2
%                          Jelin JOHN
%                 Electron content of ionosphere 

clc
clear all
close all
%%Data input

f1=1575.42e+6;
c=299792458;
lambda1=c/f1;

f2=1227.60e+6;
lambda2=c/f2;
%Define the position of the observer
x10=[4.159404458308991e+6,672972.065,4.77245255894603e+6];
x20=[4.1551687043802983e+6,672949.963,4.776113966971923e+6];

%Input phase data,extract necessay data and Reshape data according to different satellites
phase1=load('Phase1.dat');
phase2=load('Phase2.dat');
% phase1=phase1(:,1);
% phase2=phase2(:,1);
% phase1=reshape(phase1,11,25);
% phase2=reshape(phase2,11,25);

%Input visible satellite position
visSat=load('VisibleSatellites.dat');
% visSat_x=reshape(visSat(:,1),11,25);
% visSat_y=reshape(visSat(:,2),11,25);
% visSat_z=reshape(visSat(:,3),11,25);

%%
%Single differences
delta_phase=phase1(:,1)-phase2(:,1);
t=phase1(:,2);
k=sum(t==0);%To aquire the number of satellites
epoch=length(phase1(:,1))/k;

%Double differences]
D=[eye(k-1,k-1),-ones(k-1,1)];
delta_phi=D*delta_phase(1:k);
for i=1:epoch-1
   delta_phii=D*delta_phase((i*k+1):(i*k+k));
   delta_phi=[delta_phi;delta_phii];
end

P=(1/k)*(-ones(k-1,k-1)+diag(k*ones(k-1,1)));
P_bar=P;
for i=1:epoch-1
    P_bar=blkdiag(P_bar,P);
end

%Compute the measurement range by phase measurement
[M,N]=size(visSat(:,1:3));
relvec1=visSat(:,1:3)-repmat(x10,M,1);
relvec2=visSat(:,1:3)-repmat(x20,M,1);

phase10=zeros(length(visSat),1);
for i=1:length(visSat)
   phase10(i)=sqrt(relvec1(i,1:3)*relvec1(i,1:3)');
end

phase20=zeros(length(visSat),1);
for i=1:length(visSat)
    phase20(i)=sqrt(relvec2(i,1:3)*relvec2(i,1:3)');
end

delta_phase0=phase10-phase20;
delta_phi0=D*delta_phase0(1:k);
for i=1:epoch-1
   delta_phiii=D*delta_phase0((i*k+1):(i*k+k));
   delta_phi0=[delta_phi0;delta_phiii];
end


Y=[delta_phi-delta_phi0];

nrm=sqrt(diag(relvec1(:,1:3)*relvec1(:,1:3)'));
el1=diag(1./nrm)*relvec1(:,1:3);

el_diff1=D*el1(1:k,:);
for i=1:epoch-1
    el_diffi=D*el1((i*k+1):(i*k+k),:);
    el_diff1=[el_diff1;el_diffi];
end
el_diff1=-el_diff1;

lam1=lambda1*eye(k-1,k-1);
for i=1:epoch-1
lam10=lambda1*eye(k-1,k-1);
lam1=[lam1;lam10];
end

A=[el_diff1,lam1];
beta_hat= inv(A'*P_bar*A)*A'*P_bar*Y

%%
%Approximate the matrix Q22 by its main diagonal
beta_hat1=beta_hat(1:3);
beta_hat2=beta_hat(4:3+k-1);

A1=A(:,1:3);
A2=A(:,4:3+k-1);

%Looking for N and then Q
N=[A1'*P_bar*A1,A1'*P_bar*A2;A2'*P_bar*A1,A2'*P_bar*A2];
Q=inv(N);
Q22_ori=Q(4:(3+k-1),4:(3+k-1));
Q22=diag(Q22_ori)
Q22=diag(Q22);
[v,w]=eig(Q22);
max_w=max(diag(w));

% 

%Compute the confidence region
sigma_square=(Y-A*beta_hat)'*P_bar*(Y-A*beta_hat)/(250-13)
sigma=sqrt(sigma_square);
lon_axi=sqrt(sigma_square*max_w*284.66)

beta2=round(beta_hat2)

%Check
kapa=284.66;
kapacheck=(beta2-beta_hat2)'*inv(Q22)*(beta2-beta_hat2)
kapa*sigma_square
%Resolve the ambiguity double differences
beta_hat_1i=inv(A1'*P_bar*A1)*A1'*P_bar*(Y-A2*beta2)




