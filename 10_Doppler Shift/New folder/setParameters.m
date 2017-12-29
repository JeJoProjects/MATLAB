%%%%%%  sets parameters

global f;
global L;
global If;
global fs;
global dt;
global df;

f=1.57542e+9;   % L1 frequency
L=1.0e-3/1023;  % chip length
If=1.023e+6;    % intermediate frequency
fs=4*If;        %  sampling frequency
dt=1.0/fs;      %  sampling width
df=500;         %  Doppler step

