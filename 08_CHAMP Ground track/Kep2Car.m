function [ X ] = Kep2Car( elements, epochs )

GM=3.986005e+14;
n=sqrt(GM/elements(1)^3);
M=n*epochs;

%%% iterations

E0=M;
for i=1:6
    E1=E0+ ( elements(2)*sin(E0) );
    E0=E1;    
end

XietaZeta=elements(1)*[ cos(E1)-elements(2);...
                       sqrt(1-elements(2)^2)*sin(E0);...
                       zeros(size(E1))];
                   
x=R3( -elements(4)*R1( -elements(3)*R3(-elements(5)) ) )*XietaZeta;
end

