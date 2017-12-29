clc
clear all
close all

a=6377397.155;
f1=299.1528153513233;
ex= (2*f1 )- (f1^2);

[Station, B, L, H]= textread('PotsdamLBH.dat','%s %f %f %f')
[Station, X, Y, Z]= textread('GREFXYZ.dat','%s %f %f %f')

%%% Conversion L,B,H-->X,Y,Z
[x,y,z]=LBH2xyz(B*pi/180,L*pi/180,H);
delta=[x-X,y-Y,z-Z];
format bank
checkk=sqrt(diag(delta'*delta));

A1=repmat(eye(3),3,1)
aux=[X';Y';Z'];
A2=aux(:)
shapeMatrix=[0 0 0;...
             0 0 1;...
             0 -1 0;...
             0 0 -1;...
             0 0 0 ;...
             1 0 0;...
             0 1 0;...
             -1 0 0;...
             0 0 0];
aux=shapeMatrix*aux;
A3=reshape(aux,3,[])';
A=[A1, A2, A3]
Yv=Delta'
beta=A\Yv(:)




