function [ Uprime ] = ODE( t,u )

global c2;
global h;


u1=u(1:401);  % U prime
u2=u(402:802);% % U double prime
a=conv([1 -2 1],u1, 'full' );

up1=u2;        % U prime
up2=(c2/h^2).*a(2:end-1);
Uprime=[up1; up2];

end

