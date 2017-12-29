%%%% single difference phase

clc
clear all

close all
%%%% parameters

f1=1575.42e+6;
f2=1227.60e+6;

x10=[4.159404458308991e+6, 672972.0656136316, 4.77245255894603e+6];
x20=[4.155168740380298e+6, 672949.9632179659, 4.776113966971923e+6];

alpha1=f1^2/(f1^2-f2^2);
alpha2=-f1^2/(f1^2-f2^2);
%%%% read data

%%%%% data load

PR11=load('PR11OrderedByTime.dat');
PR12=load('PR12OrderedByTime.dat');
PR21=load('PR21OrderedByTime.dat');
PR22=load('PR22OrderedByTime.dat');

SatPos=load('VisibleSatellites.dat');

Phi1=load('PhaseRec1OrderedByTime.dat');
Phi2=load('PhaseRec2OrderedByTime.dat');


%%%%%% residual pseudoranges

relvec1=[SatPos(:,1:3)-repmat(x10,length(SatPos),1),SatPos(:,5)];
relvec2=[SatPos(:,1:3)-repmat(x20,length(SatPos),1),SatPos(:,5)];

PR10=sqrt(diag(relvec1(:,1:3)*relvec1(:,1:3)'));
PR20=sqrt(diag(relvec2(:,1:3)*relvec2(:,1:3)'));

% figure
% plot(PR11(1:11:220,1)-PR10(1:11:220))
% title('influence of code range errors for satellite 1')

DPR11=[PR11(:,1)-PR10,PR11(:,2)]; %% receiver 1 frequency 1
DPR12=[PR12(:,1)-PR10,PR12(:,2)]; %%% receiver 1 ferquency 2

DPR21=[PR21(:,1)-PR20,PR21(:,2)]; %% receiver 2 frequency 1
DPR22=[PR22(:,1)-PR20,PR22(:,2)]; %%% receiver 2 ferquency 2



%%%%%%% L3 solutions

DPRL31=[alpha1*DPR11(:,1)+alpha2*DPR12(:,1),DPR11(:,2)];

beta=SinglePos(DPRL31,relvec1,0)';   %%%% very bad style
for t=60:60:1140
    beta1=SinglePos(DPRL31,relvec1,t)';
    beta=[beta;beta1];
end

figure
plot(beta(:,4))
title('clock error of receiver 1')
dt1=mean(beta(:,4));

figure
plot3(beta(:,1),beta(:,2),beta(:,3),'.','Markersize',10)
title('L3 point scatter of receiver 1');

DPRL32=[alpha1*DPR21(:,1)+alpha2*DPR22(:,1),DPR21(:,2)];

beta=SinglePos(DPRL32,relvec2,0)';   %%%% very bad style
for t=60:60:1140
    beta1=SinglePos(DPRL32,relvec2,t)';
    beta=[beta;beta1];
end

figure
plot(beta(:,4))
title('clock error of receiver 2')
dt2=mean(beta(:,4));

figure
plot3(beta(:,1),beta(:,2),beta(:,3),'.','Markersize',10)
title('L3 point scatter of receiver 2');



%%%%% SD phase solution with estimated clocks

DeltaPhi=[Phi1(:,1)-Phi2(:,1),Phi1(:,2)];
DeltaPhiR=[DeltaPhi(:,1)-(PR10-PR20)-(dt1-dt2)*ones(length(PR10),1),Phi1(:,2)];

beta=SingleDiffPhase(DeltaPhiR,relvec1)


%%% SD Phase solution with real clocks

DeltaPhiR=[DeltaPhi(:,1)-(PR10-PR20)-370.0*ones(length(PR10),1),Phi1(:,2)];

beta0=SingleDiffPhase(DeltaPhiR,relvec1)