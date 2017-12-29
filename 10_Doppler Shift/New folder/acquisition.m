%%%%% acquisition

clc 
close all
clear all


setParameters

RecSig=load('ReceivedSignal.dat');

plot(RecSig(1:60))

t=0:dt:0.001;

figure
plot(t(1:60),CRshiftet(t(1:60),0,0),t(1:60),CRshiftet(t(1:60),L/3,0))

CA=load('CACode.dat');
figure

plot(t(1:60),CAshifted(t(1:60),0,0,CA),t(1:60),CAshifted(t(1:60),1.5*L,0,CA))


%%%%%% time ascquisiion


%y=TimeAcq(RecSig,CA);

%figure
%mesh(abs(y))

 %%%%% frequency acquisition
 
 tic
y=FreqAcq(RecSig,CA);
toc

figure
mesh(abs(y))

