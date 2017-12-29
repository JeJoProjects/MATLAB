%%%% Singel differences
clc
clear all

close all
%%%% parameters

f1=1575.42e+6;
f2=1227.60e+6;

x10=[4.159404458308991e+6, 672972.0656136316, 4.77245255894603e+6];
x20=[4.155168740380298e+6, 672949.9632179659, 4.776113966971923e+6];


%%%% read data

%%%%% data load

PR11=load('PR11OrderedByTime.dat');
PR12=load('PR12OrderedByTime.dat');
PR21=load('PR21OrderedByTime.dat');
PR22=load('PR22OrderedByTime.dat');

SatPos=load('VisibleSatellites.dat');


trueDX=load('trueDX.dat');


%%%%%% single point positioning on f1 for both receivers

dx1=zeros(3,length(PR11(:,2)));
clock1=zeros(1,length(PR11(:,2)));


for i=1:1426
    t=(i-1)*60;

    [beta, GDOP]=SPP_epoch(PR11,SatPos,t,x10);
    dx1(:,i)=beta(1:3);
    clock1(i)=beta(4);
end

figure
plot3(dx1(1,:),dx1(2,:),dx1(3,:),'.')
hold on
plot3(trueDX(1,1),trueDX(1,2),trueDX(1,3),'r.','Markersize',20)
title('sigle point for receiver 1')

dx2=zeros(3,length(PR11(:,2)));
clock2=zeros(1,length(PR11(:,2)));


for i=1:1426
    t=(i-1)*60;

    [beta, GDOP,rho0]=SPP_epoch(PR21,SatPos,t,x20);
    rho0
    dx2(:,i)=beta(1:3);
    clock2(i)=beta(4);
end

figure
plot3(dx2(1,:),dx2(2,:),dx2(3,:),'.')
hold on
plot3(trueDX(2,1),trueDX(2,2),trueDX(2,3),'r.','Markersize',20)
title('single point for receiver 2')
%%%%%%% single difference solution

ddx=zeros(3,length(PR11(:,2)));
dclock=zeros(1,length(PR11(:,2)));


for i=1:1426
    t=(i-1)*60;

    beta=SD_epoch(PR11,PR21,x10,x20,SatPos,t);
    ddx(:,i)=beta(1:3);
    dclock(i)=beta(4);
end

trueB=trueDX(1,:)-trueDX(2,:);
figure
plot3(ddx(1,:),ddx(2,:),ddx(3,:),'.')
hold on
plot3(trueB(1),trueB(2),trueB(3),'r.','Markersize',20)
title('single difference solution')
%%%%%%% L3 solutions

alpha1=f1^2/(f1^2-f2^2);
alpha2=-f2^2/(f1^2-f2^2);

PR31=[alpha1*PR11(:,1)+alpha2*PR12(:,1), PR21(:,2)];


dx31=zeros(3,length(PR11(:,2)));
clock31=zeros(1,length(PR11(:,2)));


for i=1:1426
    t=(i-1)*60;

    [beta, GDOP]=SPP_epoch(PR31,SatPos,t,x10);
    dx31(:,i)=beta(1:3);
    clock31(i)=beta(4);
end


PR32=[alpha1*PR21(:,1)+alpha2*PR22(:,1), PR12(:,2)];
dx32=zeros(3,length(PR12(:,2)));
clock32=zeros(1,length(PR12(:,2)));


for i=1:1426
    t=(i-1)*60;

    [beta, GDOP]=SPP_epoch(PR32,SatPos,t,x20);  %%%%% x20 istead of x10 !!!!
    dx32(:,i)=beta(1:3);
    clock32(i)=beta(4);
end


db=(dx31-dx32);

figure
plot3(db(1,:),db(2,:),db(3,:),'.')
hold on
plot3(trueB(1),trueB(2),trueB(3),'r.','Markersize',20)
title('difference of L3 solutions for both receivers')