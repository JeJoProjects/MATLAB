%% Ionosheric Range Correction
clc; clear all; close all;

SatPos=load('SatPos.dat'); 
X=[4.159425301753478e+6,672934.917, 4.772423756493888e+6]; %% Position vec. of OBSERVER.
PosVec_Sat=(1000*SatPos(:,2:4))-repmat(X,3072,1); %% Position vec. of ALL satallites.
% X2=[4.159404458308991e+6, 672972.065, 4.77245255894603e+6];

norm_X=X/norm(X);  %% normalized position vector of OBSERVER.
%%% normalized position vector ALL satellites.
norm_PosVec_Sat=PosVec_Sat./repmat(sqrt(sum([PosVec_Sat.*PosVec_Sat]'))',1,3);

%% Cos(angle b/w the two vectors) in a matrix for all satellites.
cosz=norm_PosVec_Sat*norm_X';

%%% To find the intex of visible satellites and save it as a verible.
VisSat=SatPos(  find(cosz(:)>cosd(80)) ,:);
% figure; plot3(VisSat(:,2),VisSat(:,3),VisSat(:,4));
% title(' All Visible Sats.');

Sat10=SatPos( find(SatPos(:,1)==10) ,:); %% To find& save only Satellite 10's Value in a variable.
VisSat10=SatPos( find(VisSat(:,1)==10) ,:); %% To find& save only Visible portions of Satellite 10 in a variable.
% figure; plot3(VisSat10(:,2),VisSat10(:,3),VisSat10(:,4),'r','linewidth',2.5);
% hold on; 
% grid on;
% plot3(Sat10(:,2),Sat10(:,3),Sat10(:,4),'k','linewidth',1.5);
% title(' Visible Portion of  Sat.= 10');

%% Calculation of TEC. for Visible Satellite 10
R=6371e+3; 
H=1000e+3; 
VTEC=86e+16;

sinz=sqrt(1-cosz.^2);
SinP=(R/(R+H))*sinz;
CosP=sqrt(1-SinP.^2); 

TEC=VTEC./CosP;






