clc
clear all
close all

%%% Set parameters
RecSig= load('ReceivedSignal.dat');
CA=load('CACode.dat');
t=0:dt:0.001;

figure; plot(RecSig(1:60));

figure; plot(t(1:60),CRshifted(t(1:60),0,0), t(1:60),CRshifted(t(1:60),L/3,0));

figure; plot(t(1:60),CAshifted(t(1:60),0,0), t(1:60),CAshifted(t(1:60),1.5*L,0));

%%% Time Acq

y=TimeAcq(RecSig,CA);
figure; plot();