function [ x,y,z] = LBH2xyz( B,L,H )

a=6377397.155;
f1=299.1528153513233;
ex= sqrt((2*f1)- (f1^2));

N=a./sqrt(1-ex^2*sin(B).^2);
x=(N+H).*cos(B).*cos(L);
y=(N+H).*cos(B).*sin(L);
z=((1-ex^2).*N+H).*cos(L);
end

