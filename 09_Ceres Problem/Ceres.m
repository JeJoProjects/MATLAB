%%%%% Ceres Problem

clc
close all

%%%%%% data input
global GM

GM=3.982e+14;

L=9.16*pi/180;
B=48.75*pi/180;


H=300;

a=6378137.0;
e2=6.943800229e-3;
GAST0=(23.0+ 9.0/60)*pi/12;  %%%%% wrong GAST0
om=270*pi/180;               %%%%% wrong om
omega=7.292115e-5;
ex=0.3;

%%%%% observer position

Nr=a/sqrt(1-e2*sin(B)^2);

x80=[(Nr+H)*cos(B)*cos(L);(Nr+H)*cos(B)*sin(L);(Nr*(1-e2)+H)*sin(B)]; %%verified


%%%%%% read relative vectors

dat=dlmread('CeresObs.dat','',1);

d=dat(:,2:4)';

%%%%% WGS84 positions

xWGS84=repmat(x80,1,20)+d;  %% verified



%%%%%% transformation to ICRS

GAST=GAST0+omega*dat(:,1);
size(GAST)
xICRS=[cos(GAST)'.*xWGS84(1,:)-sin(GAST)'.*xWGS84(2,:);...
      sin(GAST)'.*xWGS84(1,:)+cos(GAST)'.*xWGS84(2,:);...
      xWGS84(3,:)];   %%%% verified
  
  
%%%%%%%%  normal vector


nv=cross(xICRS(:,1),xICRS(:,20));  %%%% x!!!ICRS and not d 
nvn=nv/norm(nv);  %%%%% verified



inc=acos(nvn(3))*180/pi;  %%% verified

Omega=(atan2(nvn(2),nvn(1))+pi/2)*180/pi;   %%%% verified


%%%%%% transformation to perifocal system

xietazeta=R3(om)*R1(inc*pi/180)*R3(Omega*pi/180.0)*xICRS;



figure
plot(xietazeta(1,:),xietazeta(2,:))


%%%%%% semimajor axis
nu=atan2(xietazeta(2,:),xietazeta(1,:));
r=sqrt(diag(xietazeta'*xietazeta))




at=r'.*(1+ex*cos(nu))/(1-ex^2);

figure
plot(at)
aest=mean(at);

%%%% M0
 sinE0=sqrt(1-ex^2)*sin(nu(1))/(1+ex*cos(nu(1)));
 E0=asin(sinE0);
 M0=(E0-ex*sinE0)*180/pi;
 
 aest
 ex
 inc
 Omega
 om
 M0
 
 %%%%% check the solution
 
 el=[aest,ex,inc*pi/180,Omega*pi/180,om,M0*pi/180];
  t=dat(:,1)';
  
  orbICRS=Kep2Cart(el,t);
  
  err=orbICRS-xICRS;
  
  figure
  plot3(orbICRS(1,:),orbICRS(2,:),orbICRS(3,:))
  hold on
  plot3(xICRS(1,:),xICRS(2,:),xICRS(3,:),'r')