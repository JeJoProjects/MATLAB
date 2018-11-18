% Lab1: Astronomical position determination (Standlinie)
% Jelin John 

clc; clear all; close all;
%% Given that:
yy=2014; mm=11;dd=27; hh=19;

alpha1=(17+46/60+1.94/3600)*pi/12;  alpha2=(18+24/60+6.44/3600)*pi/12; 
delta1=56.73*pi/180;                delta2=37.9*pi/180;
z1=52.9823*pi/180;                  z2=59.3379*pi/180;

thetap=48.75*pi/180;  lambdap=9.16*pi/180;

%% parameters (Semi Major axis (a) & Flattening (f))
a1=6378137      ; f1=1/298.257;     % GRS80
a2=6377397.155 ; f2=1/299.153;    % Bessel1841

%% Unit vectors in ICRS for star 1 & star 2
ICRS1=[cos(delta1)*cos(alpha1);cos(delta1)*sin(alpha1);sin(delta1)];
ICRS2=[cos(delta2)*cos(alpha2);cos(delta2)*sin(alpha2);sin(delta2)];

%%  Precession correction
T=(G2JD(yy,mm,dd,hh,0,0)-2451545.0)/36525.0;
 
zetaA=(2306.2181 + 1.39656*T - 0.000139*T^2) *T...
                              + (0.30188 - 0.000344*T) *T^2 +0.017998* T^3;
zetaA=zetaA/3600*pi/180;
zA=(2306.2181 + 1.39656*T - 0.000139*T^2)* T...
                             + (1.09468 - 0.000066* T)* T^2 +0.018203* T^3;
zA=zA/3600*pi/180;
thetaA=(2004.3109 -0.85330* T - 0.000217* T^2) *T ...
                               - (0.42665-0.000217* T)* T^2 -0.041833* T^3;
thetaA=thetaA/3600*pi/180;
 
CEP1=R3(-zA)*R2(thetaA)*R3(-zetaA)*ICRS1; 
CEP2=R3(-zA)*R2(thetaA)*R3(-zetaA)*ICRS2; 
 
%% Earth rotation correction
[HH MM SS]=GMST(yy,mm,dd,hh,0,0);  
 
THETA=(HH+MM/60+SS/3600)*pi/12 ; 

%% unit vectors in ITRS 
ITRS1=R3(THETA)*CEP1;
ITRS2=R3(THETA)*CEP2;

Tta_1=atan2(ITRS1(3,1),(sqrt(ITRS1(1,1)^2+ITRS1(2,1)^2)));
Tta_2=atan2(ITRS2(3,1),(sqrt(ITRS2(1,1)^2+ITRS2(2,1)^2)));

Lda_1=atan2(ITRS1(2,1),ITRS1(1,1));
Lda_2=atan2(ITRS2(2,1),ITRS2(1,1));

%% Ellipsoid Co-ordinates
e1=2*f1-f1^2; e2=2*f2-f2^2;
N=a1/sqrt(1+e1^2*sin(thetap*pi/180));

XYZ=[N*cos(thetap)*cos(lambdap); N*cos(thetap)*sin(lambdap); N*(1+e1^2)*sin(thetap)];
rot=[1+11.99*10^-6 (3.3964/3600)*pi/180 -((0.5355/3600)*pi/180); -(3.3964/3600)*pi/180 1+11.99*10^-6 (-1.0778/3600)*pi/180; (0.5355/3600)*pi/180 -(-1.0778/3600)*pi/180 1+11.99*10^-6];
XYZ_NEW=[-588.196;-108.790;-378.506]+rot*XYZ;

%% Position of the observer with respect to Rauenberg datum
L_final=atan2(XYZ_NEW(2),XYZ_NEW(1));
B_final=(atan2(XYZ_NEW(3)*sin(L_final),(1-e2^2)*XYZ_NEW(2)))*180/pi
L_final=L_final*180/pi
