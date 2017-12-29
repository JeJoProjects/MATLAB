clc
clear all
close all

a=6679137;    %%%% Keplerian Elemands
e=0.000325;
i=87.218    *(pi/180);
Ohm=229.0417*(pi/180);
Omega=83.3546*(pi/180);
M=276.8105   *(pi/180);
GM=3.986005e+14;

format long
elements=[a, e, i, Ohm, Omega, M];

%%% Possibility Check

T=2*pi*sqrt(a^3/GM);

epochs=linspace(0,T);
x=Kep2Car(elements,epochs);

figure
plot3(x(1,:),x(2,:),x(3,:))
title('Orbits in CIRS');


tlaunch= G2JD(2000,7,15,11,59,59);
tstar  = G2JD(2009
%%%% non-real data



% x1=Kep2Car(elements,0.0);
% x2=Kep2Car(elements,T);
% closing=x2-x1;
