

clear all
clc
close all

format long

%%%% data input

el=[6679137.0,0.000325,87.218,229.0417,83.3546,276,8105];

%%%% convert to rad

factor=pi/180*ones(size(el));
factor(1)=1.0;
factor(2)=1.0;

elc=el*diag(factor);

%%%% plausibility check

GM=3.986005e+14;
T=2*pi/sqrt(GM/elc(1)^3);
epochs=linspace(0,T);
x1 = Kep2Car( elc,epochs );

figure
plot3(x1(1,:),x1(2,:),x1(3,:))
title('orbit in CIRS')


%%%% now real data

tlaunch=G2JD(2000,7,15,11,59,59);
tstart=G2JD(2000,7,1,0,0,0);
tstop=G2JD(2000,7,31,0,0,0);

epochs=linspace(tstart-tlaunch,tstop-tlaunch,7200)*86400

x = Kep2Car( elc,epochs );

figure
plot3(x(1,:),x(2,:),x(3,:))
title('orbit in ICRS')

%%%% transformation to ITRS

GAST0=GMST(2009,7,1,0,0,0);
GAST=GAST0*pi/12+2*pi*linspace(0,tstop-tstart,7200)

xCTS=[cos(GAST).*x(1,:)+sin(GAST).*x(2,:);-sin(GAST).*x(1,:)+cos(GAST).*x(2,:);x(3,:)];

figure
plot3(xCTS(1,:),xCTS(2,:),xCTS(3,:))
title('orbit in ITRS')


%%%% ground track

L=atan2(xCTS(2,:),xCTS(1,:))*180/pi;
L=L';
B=asin(xCTS(3,:)./sqrt(xCTS(1,:).^2+xCTS(2,:).^2+xCTS(3,:).^2))*180/pi;
B=B';

figure
load('coast')
plot(L,B)
hold on
plot(long,lat,'r')
title('ground track')

L1=circshift(L,1);
DL=L-L1;

pos=find(abs(DL)>200);
L(pos)=NaN;
B(pos)=NaN;

figure
plot(L,B)
hold on
plot(long,lat,'r')
title('ground track')