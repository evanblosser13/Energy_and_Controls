%% AME 4393/5393-003: Renewable Energy Systems and Control
%   HW 2 : State Estimation in Oil Well Drilling
%   Evan BLosser, Sooner ID 113489794
%   Prof. Dong Zhang
%   Due May 1, 2024

% Sooner_Boomer_HW5.m
clear; close all;
fs = 15;    % Font Size for plots

%% Drill String Parameters

J_T = 100;  % Table/top rotational inertia
J_B = 25;   % Bottom/bit rotational inertia
k = 2;     % Spring constant
b = 5;     % Drag coefficient                 

%% Problem 2 - Observability Analysis
% State space matrices
A4 = [-b/J_T, 0, -k/J_T, k/J_T; 
       0, -b/J_B, k/J_B, -k/J_B; 
       1, 0, 0, 0; 
       0, 1, 0, 0];

B4 = [1/J_T, 0;
       0, -1/J_B;
       0, 0;
       0, 0];

C4 = [1, 0, 0, 0];

% Compute observability Matrix for 4-state system and rank
O4 = obsv(A4, C4);
disp('Rank of Observability Matrix for four-state system');
disp(rank(O4));

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
O = obsv(A, C);
disp('Rank of Observability Matrix for three-state system');
disp(rank(O));

%% Problem 3 - Measurement Data
data = csvread('HW5_Data.csv');

t = data(:,1);       % t   : time vector [sec]
y_m = data(:,2);     % y_m : measured table velocity [radians/sec]
T= data(:,3);        % T: table torque [N-m]
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

% Set Gain scale (2-10)
Luen_s = 3;
% move eigen values
move =  -0.02;

% Eigenvalues of open-loop system
disp('Eigenvalues of 3-state system:')
Eg_A = eig(A)
S1=real(Eg_A(1));
S2=real(Eg_A(2));
S3=real(Eg_A(3));
if S1 < 0 && S2 < 0 && S3 < 0
    Stability = "Stable";
else
    Stability = "Marginally/Unstable";
end
disp('-----------Stability-----------')
disp("Lambda = " + S1 + ", " + S2 + ", " + S3);
disp("The system is " + Stability);
% Desired poles of estimation error system
%   They should have negative real parts
%   Complex conjugate pairs
Re_Pole   = (Eg_A(1) + move )*Luen_s;
i_pole    = (Eg_A(2) + move )*Luen_s;
i_pole_n  = (Eg_A(3) + move )*Luen_s;

% Assign eigen values
lam = [Re_Pole, i_pole,  i_pole_n];

% Compute observer gain 
% (See Remark 3.1 in Notes. Use "place" command)
L = place(A',C', lam)' 
disp(L);

% State-space Matrices for Luenberger Observer
A_lobs = (A - L*C);
B_lobs = [B - L*D, L]; 
C_lobs = C;
D_lobs = [0,0];

% Inputs to observer
u = [T, y_m];

% Initial Conditions
x_hat0 = [0,0,0];

% Simulate Response
sys_lobs = ss(A_lobs, B_lobs, C_lobs, D_lobs);
[y, tsim, x_hat] = lsim(sys_lobs, u, t, x_hat0);

% Parse states
theta_hat = x_hat(:,3);
omega_T_hat = x_hat(:,1);
omega_B_hat = x_hat(:,2);

% Calculate RMS
luen_est_err = omega_B_true - omega_B_hat;
RMSE = sqrt(mean((omega_B_true - omega_B_hat).^2));
disp(['Luenberger RMSE: ', num2str(RMSE), ' rad/s']);


% Plot Results
figure(2);
subplot(2,1,1);
% Plot true and estimated bit velocity
plot(t, omega_B_hat,'LineWidth',2)
hold on
plot(t, omega_B_true,'LineWidth',2)
hold off
xlabel('Time [sec]','FontSize',fs)
ylabel('Bit Velocity: w(t) [radians/sec]','FontSize',fs)
title('True vs. Estimated Bit Velocity')
legend('Estimated', 'True')
set(gca,'FontSize',fs);



% Calculated error
luen_est_err = omega_B_true - omega_B_hat;
RMSE = sqrt(mean((omega_B_true - omega_B_hat).^2));
disp(['Luenberger RMSE: ', num2str(RMSE), ' rad/s']);

Error_Estim = omega_B_true - omega_B_hat;

subplot(2,1,2);
% Plot error between true and estimated bit velocity
plot(t, Error_Estim,'LineWidth',2)
xlabel('Time [sec]','FontSize',fs)
ylabel('Error [%]','FontSize',fs)
set(gca,'FontSize',fs);

% % Problem 5 - Kalman Filter
% 
% % Noise Covariances
% W = 0.0005 * eye(3);
% N = 0.02;
% Sig0 = eye(3);
% 
% % Initial Condition
% x_hat0 = [0,0,0];
% states0 = [x_hat0, Sig0(:)'];
% 
% 
% 
% 
% % Integrate Kalman Filter ODEs
% [t, z] = ode45(@(t,z) ode_kf(t,z,t, Torq, y_m, A, B, C, W, N), t, states0);
% 
% % Parse States
% theta_hat = z(:,3);
% omega_T_hat = z(:,1);
% omega_B_hat = z(:,2);
% Sig33 = z(:,11);
% 
% omega_B_tilde = omega_B_true - omega_B_hat;
% omega_B_hat_upperbound = omega_B_hat + sqrt(Sig33);
% omega_B_hat_lowerbound = omega_B_hat - sqrt(Sig33);
% 
% RMSE = sqrt(mean(omega_B_tilde.^2));
% disp(['Kalman Filter RMSE: ', num2str(RMSE), ' rad/s']);
% 
% % Plot Results
% figure(3);
% subplot(2,1,1);
% % Plot true and estimated bit velocity
% plot(t, omega_B_true, 'C0', t, omega_B_hat, 'C1', t, omega_B_hat_upperbound, 'C3--', t, omega_B_hat_lowerbound, 'C3--');
% xlabel('Time [sec]');
% ylabel('Bit Velocity [rads/sec]');
% title('True vs Estimated Bit Velocity');
% legend('True Bit Velocity', 'Est. Bit Velocity', 'Upper STD Bound', 'Lower STD Bound');
% 
% subplot(2,1,2);
% % Plot error between true and estimated bit velocity
% plot(t, omega_B_tilde, 'C2');
% xlabel('Time [sec]');
% ylabel('Bit Velocity [rads/sec]');
% title('True vs Estimated Error (Kalman Filter)');
% 
% % Problem 6 - Extended Kalman Filter
% 
% % New nonlinear spring parameters
% k1 = 2;
% k2 = 0.25;
% 
% % Noise Covariances
% W = 0.005 * eye(3);
% N = 0.02;
% Sig0 = eye(3);
% 
% % Initial Condition
% x_hat0 = [0,0,0];
% states0 = [x_hat0, Sig0(:)'];
% 
% % Integrate Extended Kalman Filter ODEs
% [t, z] = ode45(@(t,z) ode_ekf(t,z,t, Torq, y_m, A, B, C, W, N, k1, k2, J_T, J_B, b), t, states0);
% 
% % Parse States
% theta_hat = z(:,1);
% omega_T_hat = z(:,2);
% omega_B_hat = z(:,3);
% Sig33 = z(:,end);
% 
% omega_B_tilde = omega_B_true - omega_B_hat;
% omega_B_hat_upperbound = omega_B_hat + sqrt(Sig33);
% omega_B_hat_lowerbound = omega_B_hat - sqrt(Sig33);
% 
% RMSE = sqrt(mean(omega_B_tilde.^2));
% disp(['Extended Kalman Filter RMSE: ', num2str(RMSE), ' rad/s']);
% 
% % Plot Results
% figure(3);
% subplot(2,1,1);
% % Plot true and estimated bit velocity
% plot(t, omega_B_true, 'C0', t, omega_B_hat, 'C1', t, omega_B_hat_upperbound, 'C3--', t, omega_B_hat_lowerbound, 'C3--');
% xlabel('Time [sec]');
% ylabel('Bit Velocity [rads/sec]');
% title('True vs Estimated Bit Velocity (EKF)');
% 
% subplot(2,1,2);
% % Plot error between true and estimated bit velocity
% plot(t, omega_B_tilde, 'C2');
% xlabel('Time [sec]');
% ylabel('Bit Velocity [rads/sec]');
% title('True vs Estimated Error (EKF)');