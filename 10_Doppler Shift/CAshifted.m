function [ y ] = CAshifted( t,tau,Deltaf,CA)

global L
global If
y=CA(floor(mod((1+Deltaf/If)*(t-tau)/L)+1023,1023))+1;
end

