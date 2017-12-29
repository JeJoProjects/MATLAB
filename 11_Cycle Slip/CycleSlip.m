%%%%%% cycle slip detection

close all
clear all
clc

%%%%% constants

lambda=0.190294;

x1=[4159425.302,672934.918,4772423.756];
x2=[4155192.430,672994.323,4776103.235];


%%%% load data

Sat=load('VisSats.dat');

Phi1=load('Phase1.dat');
Phi2=load('Phase2.dat');

%%%% compute approximate slant range

vec1=Sat(:,1:3)-repmat(x1,length(Sat),1);
vec2=Sat(:,1:3)-repmat(x2,length(Sat),1);

r1=sqrt(diag(vec1*vec1'));
r2=sqrt(diag(vec2*vec2'));

% residual phase

dPhi1=Phi1-r1;
dPhi2=Phi2-r2;


rdPhi1=reshape(dPhi1,5,20);
% figure
% hold on
% for i=1:5
%     plot(rdPhi1(i,:))
% end

%%%%%% single differences

DeltaPhi=dPhi1-dPhi2;
rDeltaPhi=reshape(DeltaPhi,5,20);


% figure
% hold on
% for i=1:5
%     plot(rDeltaPhi(i,:))
% end


%%%%%% double differences
ddif=[rDeltaPhi(1,:)-rDeltaPhi(5,:);...
    rDeltaPhi(2,:)-rDeltaPhi(5,:);...
    rDeltaPhi(3,:)-rDeltaPhi(5,:);...
    rDeltaPhi(4,:)-rDeltaPhi(5,:)];
% figure
% hold on
% for i=1:4
%     plot(ddif(i,:))
% end

%triple differences

trip=diff(ddif')';


for i=1:4
    figure
    plot(trip(i,:))
end

%%%%% find cycle slip


thres=3*sqrt(8)*0.002;

col=find(abs(trip(2,:))>thres)

height=trip(2,col)/lambda
height=round(height)
    