clc;
clear all;
yy=2014;
mm=11;
dd=27;
hh=19;
 
alpha1=(18+19/60+47.7/3600)*pi/12;
delta1=56.73*pi/180;
 
alpha2=(18+57/60+52.2/3600)*pi/12;
delta2=37.9*pi/180;
 
z1=8.05554*pi/180;
z2=12.1871*pi/180;
 
[HH,MM,SS]=GMST(yy,mm,dd,hh,0,0);
 
SiderialTime=([HH+MM/60+SS/3600]*pi/12)';
 
lambda1=alpha1-SiderialTime;
lambda2=alpha2-SiderialTime
 
t=pi/180;
 
[theta,lambda]=solve(sin(t)*sin(delta1)+cos(t)*cos(delta1)*cos(1-lambda1)-cos(z1),sin(t)*sin(delta2)+cos(t)*cos(delta2)*cos(1-lambda2)-cos(z2),t)
