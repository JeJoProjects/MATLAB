%  Lab_3 (Phase Double differences Solution)##

close all; clear all; clc;
%% Data input

f1=1575.42e+6; f2=1227.60e+6;   c=299792458;
lambda1=c/f1;  lambda2=c/f2;
 %% Approximate position of the two reciver stations (Observer) 
 
x1_0=[4.159404458308991e+6,672972.065,4.77245255894603e+6];
x2_0=[4.1551687043802983e+6,672949.963,4.776113966971923e+6];

phase1=load('Phase1.dat');                % Load phases 
phase2=load('Phase2.dat');
Sat=load('VisibleSatellites.dat');        % Satellite Positions 
k=11;                                      
Number_of_Satellites=k                    % The Number of Satellites
Epochs=length(phase1(:,1))/k              % The Number of Epochs

Dir_Vec_1=Sat(:,1:3)-repmat(x1_0,275,1);          %% Directional Vectors of...
Dir_Vec_2=Sat(:,1:3)-repmat(x2_0,275,1);           % two reciver stations

% Weight Matrix = Inverse of Covariance Matrix
P=(-ones(10,10)+ diag(k*ones(k-1,1)))*(1/k); 
P_bar=P;

for i=1:Epochs-1
    P_bar=blkdiag(P_bar,P);                   % Multi epoch Weight Matrix (Block diag matrix)
end

%% (1) - Double Differences 

t=phase1(:,2);
Temp1=diag([1 1 1 1 1 1 1 1 1 1 ]);           % Creating "D Matrix" (10x11)
D=[Temp1 -ones(10,1)];                    

% Computation of Double Differences  as Matrix product of Single Differences

DELTA_PHASE=phase1(:,1)-phase2(:,1);          % Phase Single differences

DEL_PHI_1=D*DELTA_PHASE(1:k);             
for i=1:Epochs-1                           
   del_phi_2=D*DELTA_PHASE( ((i*k)+1):((i*k)+k) );
   DEL_PHI_1=[DEL_PHI_1;del_phi_2];            % Delta Phi - 1 for Y vector
end

phase_10=zeros(275,1);
for i=1:275
   phase_10(i)=sqrt(Dir_Vec_1(i,1:3)*Dir_Vec_1(i,1:3)');
end

phase_20=zeros(length(Sat),1);
for i=1:length(Sat)
    phase_20(i)=sqrt(Dir_Vec_2(i,1:3)*Dir_Vec_2(i,1:3)');
end

Delta_Phase_0=phase_10-phase_20;               % Single Differences
DELTA_PHI_0=D*Delta_Phase_0(1:k);
for i=1:Epochs-1
   Delta_Phi_2=D*Delta_Phase_0((i*k+1):(i*k+k));
   DELTA_PHI_0=[DELTA_PHI_0;Delta_Phi_2];      % Delta Phi - 0 for Y vector
end                                                              

Y=[DEL_PHI_1-DELTA_PHI_0];                     % Y - Observation Vector

nrm=sqrt(diag(Dir_Vec_1(:,1:3)*Dir_Vec_1(:,1:3)'));
el1=diag(1./nrm)*Dir_Vec_1(:,1:3);

El_DIFF_1=D*el1(1:k,:);                        % Finding Unit Directional Vector

for i=1:Epochs-1
    el_diffi=D*el1((i*k+1):(i*k+k),:);
    El_DIFF_1=[El_DIFF_1;el_diffi];
end

El_DIFF_1=-El_DIFF_1;

LAMDA_1=lambda1*eye(k-1,k-1);             
for i=1:Epochs-1
Lam_10=lambda1*eye(k-1,k-1);
LAMDA_1=[LAMDA_1;Lam_10];                       % Lamda Part of Design Matrix - A
end

A=[El_DIFF_1,LAMDA_1];                          % Design Matrix
Beta_Hat= inv(A'*P_bar*A)*A'*P_bar*Y;           % Least Square Estimates (BETA_HAT)

%% (2) - Approximate the matrix Q22 by its main diagonal

Beta_Hat_1=Beta_Hat(1:3)        % Estimates for position differences (Delta-b)
Beta_Hat_2=Beta_Hat(4:3+(k-1))  % Float estimates of the Int. Ambiguities
A1=A(:,1:3);
A2=A(:,4:3+(k-1));

N=[A1'*P_bar*A1,A1'*P_bar*A2;A2'*P_bar*A1,A2'*P_bar*A2];  % Normal Equation Matrix 
Q=inv(N);                       % Inverse of normal Euation Matrix
Q22_ori=Q(4:(3+(k-1)),4:(3+(k-1)));
Q22=diag(Q22_ori)
Q22=diag(Q22);
[v,w]=eig(Q22);
W_Max=max(diag(w));         
%% (3) - Compute the confidence region

Sigma_Square=(Y-A*Beta_Hat)'*P_bar*(Y-A*Beta_Hat)/(250-13)
Sigma=sqrt(Sigma_Square)
Long_Axis=sqrt(Sigma_Square*W_Max*284.66)
Beta_2=round(Beta_Hat_2)

Kapa=284.66;                                              % Checks                                       
Kapa_Check=(Beta_2-Beta_Hat_2)'*inv(Q22)*(Beta_2-Beta_Hat_2)
Check_2=Kapa*Sigma_Square

%% (4) - Resolve the ambiguity double differences

Beta_Hat_1_Resolved=inv(A1'*P_bar*A1)*A1'*P_bar*(Y-A2*Beta_2)
Diff =Beta_Hat_1-Beta_Hat_1_Resolved

                       %% *********************** %%




