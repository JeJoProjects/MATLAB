%%%%%%  determination of flattening

clc
close all
clear all

global GM
global R

GM=3.986005e+14;
R=6378137.0;


%%%%% data input

pos=dlmread('positions.dat','',1)';
vel=dlmread('velocities.dat','',1)';

figure
plot3(pos(2,:),pos(3,:),pos(4,:))

el=Cart2Kep(pos(2:4,:),vel(2:4,:),pos(1,:));

figure
plot(el(4,:))
title('Omega')

figure
plot(el(3,:))
title('inclination')


%%%%  determination of Omegadot

p=polyfit(pos(1,:)*86400,el(4,:),1)

figure
plot(pos(1,:),el(4,:),pos(1,:),polyval(p,pos(1,:)*86400))


Omegadot=p(1);

%%%%  determination of C20

inc=mean(el(3,:));
a=R+400000;

n=sqrt(GM/a^3);

C20=2/3*n/cos(inc)*(a/R)^2