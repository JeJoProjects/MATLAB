function y = FreqAcq(recsig,CA )

global L
global df
global dt

y=zeros(41,4093);
s=fft(recsig);
t=0:dt:0.001;
for j=-20:20
    ss=SigShifted(t,0,j*df,CA);
    ss1=fft(ss);
    y(j+21,:)=fliplr(ifft(s'.*ss1));
end

