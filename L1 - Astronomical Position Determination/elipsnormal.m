function [ne,dnedB,dNedL]=elipsnormal(L,B)

global a
global f

e2=2*f;

Ne=a/sqrt(1-e2*sin(B)^2);


sq=sqrt(a^2*(-2+2*e2-e2^2+e2*(e2-2)*cos(2*B))/(2*e2*sin(B)^2-2));

ne=[Ne*cos(B)*cos(L);Ne*cos(B)*sin(L);Ne*(1-e2)*sin(B)]/sq;
dnedB=ne;
dNedL=ne;