function y = TimeAcq(recsig,CA)

global dt
global df
global L

t=0:dt:0.001;

y=zeros(41,1023);

for j=-20:20  %%%%% doppler search
    for i=0:1022   %%%% time delay search
        ss=SigShifted(t,i*L,j*df,CA);
        y(j+21,i+1)=ss*recsig;
    end
end

end

