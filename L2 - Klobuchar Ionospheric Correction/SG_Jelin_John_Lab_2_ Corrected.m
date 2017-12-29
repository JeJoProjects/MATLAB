%% JELIN JOHN  -  2924993
%%% Assignment 2 - Satellite Geodesy

clc; clear all; format long
 
%% Read the data files
 load 'vis_orbit_z.txt';
 fid = fopen('vis_orbit_z.txt','r');
 [epoch,x,y,z,] = textread('vis_orbit_z.txt','%f %f %f %f ','headerlines',1);
 format long;
 fclose(fid);

  
%% Conversion of observer position to x,y,z
b=6356752.313;
a=6378137;
phi=48.75*pi/180;
lambda=9.16*pi/180;
N=6390424.616;
h=1.0e+6;;

X0=(N+h)*cos(phi)*cos(lambda);
Y0=(N+h)*cos(phi)*sin(lambda);
Z0=((b^2/a^2)*N+h)*sin(phi);

r= sqrt(X0^2 + Y0^2 +Z0^2);

%% Computng satellite coordinates due to rotation effects

%%%%%% The satellite coordinates are already in WGS84. No transformation to
%%%%%% the Earth-fixed system necessary


t = epoch
t = t(1:9:81)
omega = 7.292115*10e-5;


R = zeros(3,3);

xyzt = zeros(3,4)
for i = 1:9
Ri = [cos(omega*t(i))  sin(omega*t(i)) 0; ...
        sin(omega*t(i)) -cos(omega*t(i)) 0; ...
                0 0 1];
xyzM = [x(9*(i-1)+1) x(9*(i-1)+2) x(9*(i-1)+3) x(9*(i-1)+4); ...
        y(9*(i-1)+1) y(9*(i-1)+2) y(9*(i-1)+3) y(9*(i-1)+4); ...
        z(9*(i-1)+1) z(9*(i-1)+2) z(9*(i-1)+3) z(9*(i-1)+4)];
xyzi = Ri*xyzM;
xyzt = [xyzt;xyzi];
end
xyzt = xyzt(4:end,:);

Pos_Rcv=(X0*ones(9,4)),(Y0*ones(9,4)),(Z0*ones(9,4));   %%%%  invalid MATLAB satement
Pos_SV=(xyzt(1:3:27,:)),(xyzt(2:3:27,:)),(xyzt(3:3:27,:));  %%%% invalid MATLAB statement

%%%% Length from satellites to the observer
S_RS = sqrt((xyzt(1:3:27,:)- X0*ones(9,4)).^2 + (xyzt(2:3:27,:)-Y0*ones(9,4)).^2 +(xyzt(3:3:27,:)-Z0*ones(9,4)).^2);
%%%% Length from satellites to the center of the earth
S_CS = sqrt((xyzt(1:3:27,:)- zeros(9,4)).^2 + (xyzt(2:3:27,:)-zeros(9,4)).^2 +(xyzt(3:3:27,:)-zeros(9,4)).^2);
%%%% Length from the observer to the center of the earth
S_OS = r*ones(9,4);

%% Calculate the angle

%d =(S_RS.^2 + S_OS.^2 - S_CS.^2)./(2.*S_RS.*S_OS)

Angle_RS_RC = acos((S_RS.^2 + S_OS.^2 - S_CS.^2)./(2.*S_RS.*S_OS));

z  = 180 - (Angle_RS_RC*180/pi);
 
zI = asin((r./(r + 1.0e6)).*sin(z));
    
azi=z*180/pi;
el=90-(Angle_RS_RC*180/pi);
a=z;

%%--------------Klobuchar range-corrections---------------------------------

psi = (0.0137 / (el(9,4) + 0.11)) - 0.022;

phi_1 = phi + (psi * cos(a));

lambda_1 = lambda + ((psi * sin(a)) / cos(phi));

phi_m = phi_1 + (0.064 * cos(lambda_1 - 1.617));

yy=2014;
mm=9;
dd=1;
hh_1=8;
hh_2=10;

tGPS = G2JD(yy,mm,dd,hh_1,0,0) - G2JD(yy,mm,dd,hh_2,0,0);

t_IPP = (43200 * lambda_1) + tGPS;

alpha_1 = 2.6534D-08;  
alpha_2 = 2.2772D-09; 
alpha_3 = -3.5174D-07;  
alpha_4 = 5.1246D-07;

beta_1 = 1.4918D+05;  
beta_2 = 8.4820D+04; 
beta_3 = -1.5726D+06;  
beta_4 = 4.0023D+06;

s_1 = (alpha_1 * phi_m) + (alpha_2 * phi_m) + (alpha_3 * phi_m) + (alpha_4 * phi_m);

s_2 = (beta_1 * phi_m) + (beta_2 * phi_m) + (beta_3 * phi_m) + (beta_4 * phi_m);

A_1 = max(0,s_1);

P_1 = max(72000,s_2);

X_1 = (2 * 3.14 * (t_IPP - 50400)) / P_1(9,4);

F = 1.0 + (16 * ((0.53 - el).^3));

if X_1 <= +1.57
    dT_L1 = ((5 .* ((10).^(-9))) + (s_1 .* (1 - (((X_1).^2)/2) + (((X_1).^4)/24)))) .* F;
else
    if X_1 <= -1.57
        dT_L2 = ((5 .* ((10).^(-9))) + (s_1 .* (1 - (((X_1).^2)/2) + (((X_1).^4)/24)))) .* F;
    else 
        dT_L3 = (5 * ((10).^(-9)))* F;
    end
end


%%--------------VTEC---------------------------------
epoch =epoch(1:9:81);

f1 = 1575.42*1e+6;

f2 = 1227.60*1e+6;

T = 1/(40.3*(1/f1^2 - 1/f2^2));

%% P codes
S1 = P_1(1:81:9);
S2 = P_1(2:81:9);
S3 = P_1(3:81:9);
S4 = P_1(4:81:9);
S5 = P_1(5:81:9);
S6 = P_1(6:81:9);
S7 = P_1(7:81:9);
S8 = P_1(8:81:9);
S9 = P_1(9:81:9);


%%TEC computation
TEC_S1 = T*S1;
TEC_S2 = T*S2;
TEC_S3 = T*S3;
TEC_S4 = T*S4;
TEC_S5 = T*S5;
TEC_S6 = T*S6;
TEC_S7 = T*S7;
TEC_S8 = T*S8;
TEC_S9 = T*S9;

%%VTEC computation

VTEC_S1 = TEC_S1.*cos(zI(1:81:9))
VTEC_S2 = TEC_S2.*cos(zI(2:81:9))
VTEC_S3 = TEC_S3.*cos(zI(3:81:9))
VTEC_S4 = TEC_S4.*cos(zI(4:81:9))
VTEC_S5 = TEC_S4.*cos(zI(5:81:9))
VTEC_S6 = TEC_S4.*cos(zI(6:81:9))
VTEC_S7 = TEC_S4.*cos(zI(7:81:9))
VTEC_S8 = TEC_S4.*cos(zI(8:81:9))
VTEC_S9 = TEC_S4.*cos(zI(9:81:9))





%%Plots
figure;
plot(epoch,S1,epoch,S2,epoch,S3,epoch,S4,epoch,S5,epoch,S6,epoch,S7,epoch,S8,epoch,S9);
xlabel('Epoch','Fontsize',14);
ylabel('P1','Fontsize',14);
title('P code','Fontsize',16);
legend('Sat 1','Sat 2','Sat 3','Sat 4','Sat 5','Sat 6','Sat 7','Sat 8','Sat 9');

figure;
plot(epoch,TEC_S1,epoch,TEC_S2,epoch,TEC_S3,epoch,TEC_S4,epoch,TEC_S5,epoch,TEC_S6,epoch,TEC_S7,epoch,TEC_S8,epoch,TEC_S9);
xlabel('Epoch','Fontsize',14);
ylabel('TEC','Fontsize',14);
title('TEC','Fontsize',16);
legend('Sat 1','Sat 2','Sat 3','Sat 4');

figure;
plot(epoch,VTEC_S1,epoch,VTEC_S2,epoch,VTEC_S3,epoch,VTEC_S4,epoch,VTEC_S5,epoch,VTEC_S6,epoch,VTEC_S7,epoch,VTEC_S8,epoch,VTEC_S9,epoch,VTEC,epoch,VTEC_average);
xlabel('Epoch','Fontsize',14);
ylabel('VTEC','Fontsize',14);
title('VTEC','Fontsize',16);
legend('Sat 1','Sat 2','Sat 3','Sat 4','VTEC','V-a');




for i=1:size(azi,1), 
    svx(i)=el(i)*cos(a(i)); 
    svy(i)=el(i)*sin(a(i));                     %Calculate polar co-ordinates
end

polarhg([30 60])                                %Prerequisite script used to format axis
hold on
plot( svx,svy,'.r','markers',20);               %Plot satellite location
hold off
 
%% Format output 
for i=1:el,
    text(svx(i)+7,svy(i),num2str(prn(i)), 'FontSize' ,10) ; %Add PRN labels to each point
end


axis('square')
grid on;  
set(gcf, 'Color', 'w');                        %Change background of figure from grey to white
ti = get(gca,'TightInset')   ;                 %Remove extra spacing around figure
set(gca, 'LooseInset', [0,0,0,0.01]);          %Depending on the figure, you may need to add extra spacing [left bottom width height])
print( '-dtiff',  ['skyPlot'], '-r600');       %Change "-r600" to the required DPI



