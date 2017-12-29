function [ y ] = TimeAcq( RecSig, CA )

global dt;
global If;
global L;
t=0:dt:0.001;
y=zeros(41,1023);
for j=-20:20 %%% Doppler search
    for i=0:1022 %%% Time delay search
        ss=SigShifted(t,i*L,j*df,CA);
        y(j+21,i+1)=ss*RecSig;
    end
end      


end

