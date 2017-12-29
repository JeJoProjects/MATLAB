%%%%%%% ionospheric corrections

clc
close all
clear all

%%% load satellite orbits
orbits=load('SatPos.dat');

%%%% vizualization of orbit

idx17=find(orbits(:,1)==17);

orb17=orbits(idx17,:);

plot3(orb17(:,2),orb17(:,3),orb17(:,4));
title('full orbit');


%%%% visibility

obs1=[4.159404458319e+6;672972.065917;4.77245255895e+6];
relvec=1000*orbits(:,2:4)-repmat(obs1',[96*32,1])
norm_relvec=sqrt(diag(relvec*relvec'));
relvec_n=diag(1./norm_relvec)*relvec;

norm_obs=norm(obs1);

obs1_n=obs1/norm_obs;




cosz=relvec_n*obs1_n;


id_vis=find(cosz>cos(80*pi/180));

vis_sat=orbits(id_vis,:);



idx17=find(vis_sat(:,1)==17);

orb17=orbits(idx17,:);
figure
plot3(orb17(:,2),orb17(:,3),orb17(:,4))
title('visible orbit')

cosz_vis=cosz(id_vis);
plot(cosz_vis);
return;

%%%%%% TEC
VTEC=86.0e+16;
R=6378137.0;
H=1.0e+6;
sinz=sqrt(1.0-cosz_vis.^2);

sinzP=(R/(R+H))*sinz;
coszP=sqrt(1.0-sinzP.^2);

TEC=VTEC./coszP;

tag=1:length(idx17);

figure
plot(tag,VTEC*ones(1,length(idx17)),tag,TEC(idx17))


% ionospheric corrections

f1=1.57542e+9;

dI=-40.28*TEC/f1^2;

figure
plot(dI(idx17))

