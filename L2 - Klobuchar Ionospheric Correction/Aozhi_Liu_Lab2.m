clc;
clear all;
close all;

format long

%read the data
orbit=load('orbit.dat');

%define the position of the observer
% L=48.79;
% B=9.19;

L=9.19;
B=48.79;


%Convert the observer's location into catesian coordinate
[xOB yOB zOB]=LBH2XYZ(L,B,0);
xO=[xOB yOB zOB];
xOU=[xOB/sqrt(xOB^2+yOB^2+zOB^2) yOB/sqrt(xOB^2+yOB^2+zOB^2) zOB/sqrt(xOB^2+yOB^2+zOB^2)];

%Create unit vectors of the satellites with respect to the observer
Temp1=[xOB yOB zOB];
[row1 col1]=size(orbit);
Temp2=ones(row1,1);
Temp3=Temp2*Temp1;
SatVec=orbit(:,2:4)-Temp3;
SUx=SatVec(:,1)./sqrt(SatVec(:,1).^2+SatVec(:,2).^2+SatVec(:,2).^2);
SUy=SatVec(:,2)./sqrt(SatVec(:,1).^2+SatVec(:,2).^2+SatVec(:,2).^2);
SUz=SatVec(:,3)./sqrt(SatVec(:,1).^2+SatVec(:,2).^2+SatVec(:,2).^2);
SatUni=[SUx SUy SUz];

%Find observations with zenith angle<75deg
xSEU=R2(pi/2-B*pi/180)*R3(L*pi/180)*SatVec';%transform it into local level system
xZEN=xSEU(3,:);
xZENNorm=sqrt(xSEU(1,:).^2+xSEU(2,:).^2+xSEU(3,:).^2);
zenithRad=acos(xZEN./xZENNorm);
zenithDeg=zenithRad*180/pi;
index1=find(zenithDeg<90.0);%index of the satellites visible
inSatVis=ceil(index1/9.0);%PRN of the satellites visible with redundency
inSatVisNum=unique(inSatVis);%extract PRN of the satellites visible
inSatVisTab=[inSatVisNum',histc(inSatVis,inSatVisNum)']
index2=find(histc(inSatVis,inSatVisNum)==9);%index of the satellites visible in every epoch
inSatFullVis=inSatVisNum(index2)

for i=1:length(inSatFullVis)
orbVisNum(:,i)=find(orbit(:,1)==inSatFullVis(i));
 end
 orbVisNum

%Make a skyplot of visible satellites
eleDeg=90-zenithDeg;%Get elevation angle
eleRad=eleDeg*pi/180;
azimuthRad=atan2(xSEU(2,:),-xSEU(1,:));
azimuthDeg=azimuthRad*180/pi;

%Full sketch of all visible epochs
eleFull=eleDeg(index1);
AzFull=azimuthDeg(index1);
figure(1)
for i=1:1:length(inSatVis)-1
skyplot(inSatVis(i),AzFull(i),eleFull(i));
hold on
end


%satellites that are full visible
figure(2)
ele2=eleDeg(orbVisNum(:,1));
Az2=azimuthDeg(orbVisNum(:,1));
h2=skyplot2(2,Az2,ele2);
hold on

ele11=eleDeg(orbVisNum(:,2));
Az11=azimuthDeg(orbVisNum(:,2));
h11=skyplot11(11,Az11,ele11);
hold on

ele13=eleDeg(orbVisNum(:,3));
Az13=azimuthDeg(orbVisNum(:,3));
h13=skyplot13(13,Az13,ele13);
hold on

ele24=eleDeg(orbVisNum(:,4));
Az24=azimuthDeg(orbVisNum(:,4));
h24=skyplot24(24,Az24,ele24);
hold on

ele28=eleDeg(orbVisNum(:,5));
Az28=azimuthDeg(orbVisNum(:,5));
h28=skyplot28(28,Az28,ele28);
hold on

ele30=eleDeg(orbVisNum(:,6));
Az30=azimuthDeg(orbVisNum(:,6));
h30=skyplot30(30,Az30,ele30);
hold on

legend([h2,h11,h13,h24,h28,h30],'SV2','SV11','SV13','SV24','SV28','SV30');
legend('boxoff');

%%
%Klobuchar range correction
TrueAzimuth=90-azimuthDeg;
azimuth=reshape(TrueAzimuth,9,31);
elevation=reshape(eleDeg,9,31);
index3=elevation>15;
tt=1:1:9;
ttlabel=8.00:0.25:10.00;
t0=mod(239,7)*24*3600+8*3600;
t1=mod(239,7)*24*3600+10*3600;
tGPS=linspace(t0,t1,9)'*ones(1,31);
dtL1=Klobuchar(L*pi/180,B*pi/180,azimuth*pi/180,elevation*pi/180,tGPS);
dtL1Vis=dtL1;
dtL1Vis(~index3)=NaN;

figure(3)
plot(ttlabel,dtL1(tt,:),'--^');grid on
xlabel('Time at UTC 2014 Sep.1st [h]');
ylabel('Klobuchar time correction for all satellies[s]');


figure(4)
plot(ttlabel,dtL1Vis(tt,:),'--^');grid on
xlabel('Time at UTC 2014 Sep.1st [h]');
ylabel('Klobuchar time correction for visible satellies[s]');
elevationVis=elevation;
elevationVis(~index3)=NaN;

figure(5)
plot(elevationVis(tt,:),dtL1Vis(tt,:),'--^');grid on
xlabel('elevation angle [deg]');
ylabel('Klobuchar time correction[s]');
