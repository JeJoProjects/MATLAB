%%%%% stand-line method
clc
clear all;

global a
global f

%%%  data input
yy=2014; mm=11; dd=27; hh=20;

alpha1=(18+19/60+47.7/3600)*pi/12;
delta1=56.73*pi/180;
alpha2=(18+57/60+52.2/3600)*pi/12;
delta2=37.9*pi/180;

z1=8.05554*pi/180;
z2=12.1871*pi/180;

a=6378137;
f=1/298.257;

L0=9.16*pi/180;
B0=48.75*pi/180;

SiderialTime=([hh+0/60+0/3600]*pi/12)';
lambda1=alpha1-SiderialTime;
lambda2=alpha2-SiderialTime;

%%%%  e cartesian in ICRS

e1ICRS=[cos(delta1)*cos(alpha1);cos(delta1)*sin(alpha1);sin(delta1)];
e2ICRS=[cos(delta2)*cos(alpha2);cos(delta2)*sin(alpha2);sin(delta2)];

%%%  transformation to ITRS
T=(juliandate(yy,mm,dd,hh,0,0)-juliandate(2000,1,0,0,0,0))/36525.0

P=Precession(T);

Theta=GMST(T)*pi/12;

e1ITRS=R3(Theta)*P*e1ICRS;
e2ITRS=R3(Theta)*P*e2ICRS;

%%%%% zenith vector

[x,y,z]=elipsnormal(L0,B0);