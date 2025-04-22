%% AME 4393/5393-003: Renewable Energy Systems and Control
%   HW 2 : State Estimation in Oil Well Drilling
%   Evan Blosser, Sooner ID 113489794
%   Prof. Dong Zhang
%   Due May 1, 2024

% Sooner_Boomer_HW5.m
clear; close all;
fs = 15;    % Font Size for plots

%% Drill String Parameters

p.J_T = 100;  % Table/top rotational inertia
p.J_B = 25;   % Bottom/bit rotational inertia
p.k = 2;      % Spring constant
p.b = 5;      % Drag coefficient

%% Problem 2 - Observability Analysis

% State space matrices
A4 = [ -p.b/p.J_T, 0, -p.k/p.J_T, p.k/p.J_T; 0, -p.b/p.J_B,p.k/p.J_B,-p.k/p.J_B; 1, 0, 0, 0; 0, 1, 0, 0];
B4 = [1/p.J_T, 0;0, -1/p.J_B;0, 0 ;0, 0];
C4 = [1, 0, 0 , 0];

% Compute observability Matrix for 4-state system and rank
O4 = obsv(A4,C4)
disp('Rank of Observability Matrix for four-state system')
rank(O4)

% New A Matrix, for 3-state system
A = [-b/J_T, 0, -k/J_T;
      0, -b/J_B, k/J_B;
      1, -1, 0];

B = [1/J_T;
      0;
      0];

C = [0, 1, 0];

D = 0; % Add empty D for Q4 

% Observability Matrix for 3-state system and rank
O = obsv(A,C)
disp('Rank of Observability Matrix for three-state system')
rank(O)

%% Problem 3 - Measurement Data
data = csvread('HW5_Data.csv');
t = data(:,1);            % t   : time vector [sec]
y_m = data(:,2);          % y_m : measured table velocity [radians/sec]
T = data(:,3);            % T   : table torque [N-m]
omega_B_true = data(:,4); % \omega_B : true rotational speed of bit [radians/sec]

figure(1); clf;

subplot(2,1,1);
% Plot table torque
plot(t, T,'LineWidth',2)
xlabel('Time [sec]','FontSize',fs)
ylabel('T(t) [N-m]','FontSize',fs)
set(gca,'FontSize',fs);

subplot(2,1,2);
% Plot measured table velocity
plot(t, y_m ,'LineWidth',2)
xlabel('Time [sec]','FontSize',fs)
ylabel('w(t) [radians/sec]','FontSize',fs)
set(gca,'FontSize',fs);


%% Problem 4 - Luenberger Observer

% Eigenvalues of open-loop system
disp('Eigenvalues of 3-state system:')
Eg_A = eig(A)

% Desired poles of estimation error system
%   They should have negative real parts
%   Complex conjugate pairs
Re_Pole   = Eg_A(1);
i_pole    = Eg_A(2);
i_pole_n  = Eg_A(3);

% Set Gain scale (2-10)
Luen_s = 5;
% Assign eigen values
lam = [Re_Pole*Luen_s, i_pole*Luen_s,  i_pole_n*Luen_s];
% Compute observer gain 
% (See Remark 3.1 in Notes. Use "place" command)
L = place(A',C', lam)' 

% State-space Matrices for Luenberger Observer
A_lobs =  A - L*C;
B_lobs =  B;
C_lobs =  C;
D_dummy = [ 0, 0;];
sys_lobs = ss( A_lobs,  B_lobs,  C_lobs,   D_dummy);

% Inputs to observer
u = [ T,  y_m];

% Initial Conditions for Luenberger Observer
x_hat0 = [ 0, 0, 0];

% Simulate Response
[y,t,x_hat] = lsim(sys_lobs,u,t,x_hat0);

% Parse out states
theta_hat   = x_hat(:,1);
omega_T_hat = x_hat(:,3);
omega_B_hat = x_hat(:,2);


% Plot Results
figure(2); clf;

subplot(2,1,1);
% Plot true and estimated bit velocity
plot(t, omega_B_hat,'LineWidth',2)
hold on
plot(t, omega_B_true,'LineWidth',2)
hold off
xlabel('Time [sec]','FontSize',fs)
ylabel('Bit Velocity: w(t) [radians/sec]','FontSize',fs)
title('True vs. Estimated Bit Velocity')
legend('omega\_B\_hat', 'omega\_B\_true')
set(gca,'FontSize',fs);


subplot(2,1,2);
% Plot error between true and estimated bit velocity
plot(t, theta_hat ,'LineWidth',2)
ylabel('Error [-]','FontSize',fs)
set(gca,'FontSize',fs);


% 
% %% Problem 5 - Kalman Filter
% % Noise Covariances
% W =  ;   % You design this one.
% N =  ;
% Sig0 = ;
% 
% % Input data
% input_data = [t, T, y_m];
% 
% % Initial Condition
% x_hat0 = [ ;  ;  ];
% states0 = [x_hat0; reshape(Sig0,[9 1])];
% 
% % Simulate Kalman Filter
% %   This coding is a little complicated for novice Matlab users.
% %   Try to reverse engineer what I've done here.
% [tsim,states] = ode45(@(t,x) ode_kf(t,x,A,B(:,1),C,input_data,W,N),t,states0);
% 
% % Parse States
% theta_hat = states(:,1);
% omega_T_hat = states(:,2);
% omega_B_hat = states(:,3);
% Sig33 =  ;   % Parse out the (3,3) element of Sigma only!
% 
% % Compute the upper and lower bounds as described in Problem 5c.
% omega_B_hat_upperbound = omega_B_hat + sqrt(Sig33); 
% omega_B_hat_lowerbound = omega_B_hat - sqrt(Sig33);
% 
% % Plot Results
% figure(3); clf;
% 
% subplot(2,1,1);
% %   Plot true and estimated bit velocity
% %   Plot estimated bit velocity plus/minus one sigma
% 
% subplot(2,1,2);
% %   Plot error between true and estimated bit velocity
% 
% 
% 
% omega_B_tilde = omega_B_true - omega_B_hat;
% RMSE = sqrt(mean(omega_B_tilde.^2));
% fprintf(1,'Kalman Filter RMSE: %1.4f rad/s\n',RMSE);
% 
% 
% %% Problem 6 - Extended Kalman Filter
% % New parameters
% p.k1 = ;
% p.k2 = ;
% 
% % now you are on your own!