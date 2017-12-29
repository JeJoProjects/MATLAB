clc; clear all; close all;

alphaMoon =(23+56/60+48.971/3600)*pi/12
deltMoon =(42.1/60+1.06/3600)*pi/180;

alphaSun=(23+58/60+1.399/3600)*pi/12;
deltSun=-(12/60+61.03/3600)*pi/180;

rSun =1.4895092174707675e+11;
rMoon=357921054;

%%%% Axix of cone

xMoon = rMoon* [cos(alphaMoon)*cos(deltMoon); sin(alphaMoon)*cos(deltMoon); sin(deltMoon)];
xSun = rSun* [cos(alphaSun)*cos(deltSun); sin(alphaSun)*cos(deltSun); sin(deltSun)];

dirVec = xMoon-xSun;

p=2*xMoon'*dirVec/(dirVec'*dirVec);
R=6378137;

q=(xMoon'*xMoon-R^2)/(dirVec'*dirVec);

t1=-p/2+sqrt(p^2/4-q);
t2=-p/2-sqrt(p^2/4-q);

xIntersect = xMoon+t2*dirVec;

norm(xIntersect)-R;


%%% rotale from celestial to terrestial by GAST
format long g;
sidAngle=GMST(2015,3,20,9,45,0);

xTerr=R3(sidAngle*pi/12)*xIntersect

%%%% Convertion to longitude and latitude

L=atan2(xTerr(2),xTerr(1))*180/pi
B=asin(xTerr(3)/R)*180/pi

