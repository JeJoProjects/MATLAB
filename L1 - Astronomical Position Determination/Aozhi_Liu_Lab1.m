clear all;
close all;
clc;
format long g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step1 xICRS

% alpha1=(18.00+19/60.00+47.7/3600.00)*pi/12.00;
% alpha2=(18.00+57/60.00+52.2/3600.00)*pi/12.00;

alpha1=(17.00+46/60.00+1.94/3600.00)*pi/12.00;
alpha2=(18.00+24/60.00+6.44/3600.00)*pi/12.00;

deta1=56.73*pi/180.00;
deta2=37.90*pi/180.00;

JD=G2JD(2014,11,27,19,0,0);

%unit vectors 
xICRS1=[cos(deta1)*cos(alpha1);cos(deta1)*sin(alpha1);sin(deta1)]
xICRS2=[cos(deta2)*cos(alpha2);cos(deta2)*sin(alpha2);sin(deta2)]

%step2 precession correction
T=(JD-G2JD(2000,1,1,12,0,0))/36525;

zetaA=(2306.2181 + 1.39656*T - 0.000139*T^2) *T...
    + (0.30188 - 0.000344*T) *T^2 +0.017998* T^3;
zA=(2306.2181 + 1.39656*T - 0.000139*T^2)* T...
    + (1.09468 - 0.000066* T)* T^2 +0.018203* T^3;
thetaA=(2004.3109 -0.85330* T - 0.000217* T^2) *T ...
    - (0.42665-0.000217* T)* T^2 -0.041833* T^3;

zetaA=zetaA/3600*pi/180;
zA=zA/3600*pi/180;
thetaA=thetaA/3600*pi/180;

xCEP1=R3(-zA)*R2(thetaA)*R3(-zetaA)*xICRS1
xCEP2=R3(-zA)*R2(thetaA)*R3(-zetaA)*xICRS2

%%%%% Earth rotation correction
[HH MM SS]=GMST(2014,11,27,19,0,0);%%GMST computation
THETA=(HH+MM/60+SS/3600)*pi/12;
xITRS1=R3(THETA)*xCEP1
xITRS2=R3(THETA)*xCEP2

%derive the position by vectors intersection
% z1=8.05554*pi/180;
% z2=12.1871*pi/180;

z1=52.9823*pi/180;
z2=59.3379*pi/180;

L0=9.16*pi/180;
B0=48.75*pi/180;

k1=cos(z1);
k2=cos(z2);

A=[xITRS1(1) xITRS1(2) xITRS1(3);xITRS2(1) xITRS2(2) xITRS2(3)];
deltaY=[k1;k2]-A*[cos(B0)*cos(L0);cos(B0)*sin(L0);sin(B0)];

B=B0;
L=L0;

C=[-sin(B)*cos(L);-sin(B)*sin(L);cos(B)];
D=[-sin(L)*cos(B);cos(B)*cos(L);0];
Anew=A*[C D];
X=Anew\deltaY;
B=B0+X(1);
L=L0+X(2);

while abs(X(2))>0.0000001
C=[-sin(B)*cos(L);-sin(B)*sin(L);cos(B)];
D=[-sin(L)*cos(B);cos(B)*cos(L);0];
Anew=A*[C D];
deltaY=[k1;k2]-A*[cos(B)*cos(L);cos(B)*sin(L);sin(B)];
X=Anew\deltaY;
B=B+X(1);
L=L+X(2);
end
L=L*180/pi
B=B*180/pi

%in cartesian coordinates

[xG yG zG]=LBH2XYZ(L,B,0);
XG=[xG yG zG]

%Transform the coordinate system into Rauenberg system
al=-1.0778/3600*pi/180;
be=0.5355/3600*pi/180;
ga=3.3964/3600*pi/180;
XR=[-588.196;-108.790;-378.506]+(1+11.99/1000000)*R1(al)*R2(be)*R3(ga)*[xG;yG;zG]
