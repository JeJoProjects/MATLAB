clc
%%%%%%%   data input

alphaMoon=(23.0+56.0/60.0+48.971/3600.0)*pi/12.0;
deltaMoon=(42.1/60.0+1.06/3600.0)*pi/180.0;
rMoon=357921054.0;

alphaSun=(23.0+58/60.0+1.399/3600.0)*pi/12.0;
deltaSun=-(12.0/60.0+61.03/3600.0)*pi/180.0;
rSun=1.4895092174707675e+11;

%%%%%% axis of shadow cone

xMoon=rMoon*[cos(alphaMoon)*cos(deltaMoon);sin(alphaMoon)*cos(deltaMoon);...
    sin(deltaMoon)];
xSun=rSun*[cos(alphaSun)*cos(deltaSun);sin(alphaSun)*cos(deltaSun);...
    sin(deltaSun)];

dirVec=xMoon-xSun


%%%%% intersection with mean Errth sphere

p= 2*xMoon'*dirVec/(dirVec'*dirVec);
R=6378137.0;

q=(xMoon'*xMoon-R^2)/(dirVec'*dirVec);

t1=-p/2.0+sqrt(p^2/4.0-q);
t2=-p/2.0-sqrt(p^2/4.0-q);

xIntersect=xMoon+t2*dirVec

%%%%%%% conversion to terectrial system

sidAngle=GMST(2015,3,20,9,45,0)

xTerr=R3(sidAngle*pi/12)*xIntersect

%%%%%% conversion to Longitude and latitude

L=atan2(xTerr(2),xTerr(1))*180.0/pi
B=asin(xTerr(3)/R)*180.0/pi

