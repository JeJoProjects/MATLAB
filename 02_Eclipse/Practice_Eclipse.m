clc
clear all
close all
alphaMoon=(23+(56/60)+(48.971/3600))*15;
deltaMoon=(42+(01/60)+(0.06/3600));
rMoon    = 357921054;
alphaSun =(23+(58/60)+(01.399/3600))*15;
deltaSun=-(12+(61/60)+(0.03/3600));
rSun    = 1.4895092174707675e+11;
cMoon=rMoon*[cosd(deltaMoon)*cosd(alphaMoon);...
             cosd(deltaMoon)*sind(alphaMoon);...
             sind(deltaMoon)];
cSun =rSun*[cosd(deltaSun)*cosd(alphaSun);...
            cosd(deltaSun)*sind(alphaSun);...
            sind(deltaSun)];        
 dirVec=cMoon-cSun


        