clc

%%%%%%  data input

alpha=(2/24.0+31/1440.0+48.7/86400.0)*2*pi
delta=(89+15/60.0+51/3600.0)*pi/180.0

L=9.16*pi/180.0;
B=48.75*pi/180.0;


%%%%% unitvector in space fixed system

e_alpha=[cos(alpha)*cos(delta);sin(alpha)*cos(delta);sin(delta)];

%%%%%  transformation to Earth fxed system

sd=GMST(2015,12,24,19,0,0)

e_terr=R3(sd*pi/12.0)*e_alpha


%%%%% transformation to horizontal system

e_hor=R2(pi/2-B)*R3(L)*e_terr

%%%%%% azimuth and elevation

el=asin(e_hor(3))*180/pi

az=atan2(e_hor(2),-e_hor(1))*180/pi