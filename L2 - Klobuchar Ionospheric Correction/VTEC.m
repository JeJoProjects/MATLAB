%%%%%%   VTEC

clear all
close all
clc

%%%%%  defining constants

f1=1575.42e+6;
f2=1227.60e+6;

L=9.19*pi/180.0;
B=48.79*pi/180.0;

R=6738137.0;

H=1.0e+6;

%%%%%  load data

orb=load('VisSats.dat');
L1=load('PR1.dat');
L2=load('PR2.dat');

%%%%%  compute TEC

TEC=1.0e-16*(L1(:,2)-L2(:,2))/(40.28*(1/f1^2-1/f2^2));

TEC=reshape(TEC,9,6);

n=1:9;

figure
plot(n,TEC(:,1),n,TEC(:,2),n,TEC(:,3),n,TEC(:,4),n,TEC(:,5),n,TEC(:,6))

title('TEC in TECU')
%%%%  computation of zenith angles

obsXYZ=R*[cos(B)*cos(L),cos(B)*sin(L),sin(B)];

obsvec=orb(:,2:4)-repmat(obsXYZ,54,1);

cosz=(obsXYZ*obsvec')./(norm(obsXYZ)*sqrt(diag(obsvec*obsvec')))';

z=acos(reshape(cosz,9,6))*180.0/pi;



figure
plot(n,z(:,1),n,z(:,2),n,z(:,3),n,z(:,4),n,z(:,5),n,z(:,6))


%%%%% converting TEC into VTEC

zprime=asin(R/(R+H)*sqrt(1.0-cosz.^2));

VTC=TEC(:).*cos(zprime)';

VTC=reshape(VTC,9,6);
figure

plot(n,VTC(:,1),n,VTC(:,2),n,VTC(:,3),n,VTC(:,4),n,VTC(:,5),n,VTC(:,6))

%%% average over satellites

mVTEC=mean(VTC');

figure
plot(n,mVTEC)
