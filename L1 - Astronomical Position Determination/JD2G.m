function [yy,mm,dd,hh,mind,sec]=JD2G(JD)
Z=floor(JD+0.5);
F=JD+0.5-Z;
if Z < 2299161 
    A=Z;
else
    S=floor((Z-1867216.25)/36524.25);
	A=Z+1+S-floor(S/4);
end
B=A+1524;
C=floor((B-122.1)/365.25);
D=floor(365.25*C);
E=floor((B-D)/30.6001);
fdd=B-D-floor(30.6001*E)+F;
if E < 14 
    mm=E-1;
else
    mm=E-13;
end
if mm >2 
    yy=C-4716; 
else
    yy=C-4715;
end
dd=floor(fdd);
fhh=(fdd-dd)*24.0;
hh=floor(fhh);
fmind=(fhh-hh)*60;
mind=floor(fmind);
sec=(fmind-mind)*60.0;