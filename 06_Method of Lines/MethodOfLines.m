%%% Methood of lines
clc
clear all
close all
global c2; 
global h;
global n;
%%% Parameters of the problem
L=2;    % which gos right and left
h=0.01; % step size
n=L/h;  % the number of lines

%%%% Initial Values
sigma=0.01; 
x=-L:h:L; size(x);

f=exp(-x.^2/sigma);
figure; 
plot(x,f); title('Initial State');grid on;

%%% Medium
 c2=[2*ones(300,1); 0.5*ones(101,1)];
 [t, u]=ode45(@ODE, 0:.01:1.2, [f'; zeros(401,1)]);
 
 figure;
 imagesc(u(:,1:401));
 
 
 figure;
 for i=1:length(t)
     plot(x,u(i,1:401))
     axis([-1.5,1.5,-0.5,1])
     grid on;
     m(i)=getframe;
 end

 
























