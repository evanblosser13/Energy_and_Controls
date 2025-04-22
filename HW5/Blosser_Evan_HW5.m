%% AME 4393/5393-003: Renewable Energy Systems and Control
%   HW 2 : State Estimation in Oil Well Drilling
%   Evan Blosser, Sooner ID 113489794
%   Prof. Dong Zhang
%   Due May 1, 2024

% Blosser_Evan_HW5.m
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

t   = data(:,1);       % t   : time vector [sec]
y_m = data(:,2);     % y_m : measured table velocity [radians/sec]
T   = data(:,3);        % T: table torque [N-m]
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
Luen_s = 4;
% move eigen values
move =  -0.24;
move_i =  -0.04;

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
Re_Pole   = (Eg_A(1) + move    )*Luen_s;
i_pole    = (Eg_A(2) + move_i  )*Luen_s;
i_pole_n  = (Eg_A(3) + move_i  )*Luen_s;

% Assign eigen values
lam = [Re_Pole, i_pole,  i_pole_n];
disp('-----Chosen eigenvalues-----')
disp((S1 + move)*Luen_s  )
disp((S2 + move_i)*Luen_s )
disp((S3 + move_i)*Luen_s  )
disp('----------------------------')
% Compute observer gain 
% (See Remark 3.1 in Notes. Use "place" command)
L = place(A',C', lam)' 

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


Error_Estim = omega_B_true - omega_B_hat;

subplot(2,1,2);
% Plot error between true and estimated bit velocity
plot(t, Error_Estim,'LineWidth',2)
xlabel('Time [sec]','FontSize',fs)
ylabel('Error [-]','FontSize',fs)
set(gca,'FontSize',fs);


%% Problem 5 - Kalman Filter

% Noise Covariances
W = 0.0042 * eye(3);
N = 0.02;
Sig0 = eye(3);

% Initial Condition
x_hat0 = [0,0,0];
states0 = [x_hat0, Sig0(:)'];


% Simulate Kalman Filter
%   This coding is a little complicated for novice Matlab users.
%   Try to reverse engineer what I've done here.
% [tsim,states] = ode45(@(t,x) ode_kf(t,x,A,B(:,1),C,input_data,W,N),t,states0);

% Integrate Kalman Filter ODEs
[tsim,states] = ode45(@(t,x) ode_kf(t,x, A, B, C,data, W, N), t, states0);

% Parse States
theta_hat   = states(:,3);
omega_T_hat = states(:,1);
omega_B_hat = states(:,2);
Sig33       = states(:,11);

Error_Estim_kf = omega_B_true - omega_B_hat;
omega_B_hat_upperbound = omega_B_hat + sqrt(Sig33);
omega_B_hat_lowerbound = omega_B_hat - sqrt(Sig33);


% Plot Results
figure(3);
subplot(2,1,1);
% Plot true and estimated bit velocity
plot(t, omega_B_true)
hold on
plot(tsim, omega_B_hat,'LineWidth',2)
hold on
plot(tsim, omega_B_hat_upperbound,'r--','LineWidth',2)
hold on
plot(tsim, omega_B_hat_lowerbound,'r--','LineWidth',2)
hold off
xlabel('Time [sec]');
ylabel('Bit Velocity [rads/sec]');
title('True vs Estimated Bit Velocity');
legend('True Bit Velocity', 'Est. Bit Velocity', 'Upper STD Bound', 'Lower STD Bound');

subplot(2,1,2);
% Plot error between true and estimated bit velocity
plot(t, Error_Estim_kf);
xlabel('Time [sec]');
ylabel('Error [-]','FontSize',fs)
title('Error (Kalman Filter)');


%% Problem 6 - Extended Kalman Filter

% New nonlinear spring parameters
k1 = 2;
k2 = 0.25;


% Noise Covariances
W = 0.005 * eye(3);
N = 0.02;
Sig0 = eye(3);


% Initial Condition
x_hat0  = [0,0,0];
states0 = [x_hat0, Sig0(:)'];

% Simulate Extended Kalman Filter  
[tsim, states_efk] = ode45(@(t,x) ode_ekf(t,x, data, k1, k2, b, J_T, J_B, N, W), t, states0);

% Parse States
theta_hat_efk   = states_efk(:,1);
omega_T_hat_efk = states_efk(:,2);
omega_B_hat_efk = states_efk(:,3);
Sig33_efk       = states_efk(:,end);

Error_Estim_efk = omega_B_true- omega_B_hat_efk;
omega_B_hat_upper_efk = omega_B_hat_efk + sqrt(Sig33_efk);
omega_B_hat_lower_efk = omega_B_hat_efk - sqrt(Sig33_efk);


% Plot Results
figure(4);
subplot(2,1,1);
% Plot true and estimated bit velocity
plot(t, omega_B_true)
hold on
plot(tsim, omega_B_hat_efk,'LineWidth',2)
hold on
plot(tsim, omega_B_hat_upper_efk,'r--','LineWidth',2)
hold on
plot(tsim, omega_B_hat_lower_efk,'r--','LineWidth',2)
hold off
xlabel('Time [sec]');
ylabel('Bit Velocity [rads/sec]');
title('True vs Estimated Bit Velocity (EKF)');
legend('True Bit Velocity', 'Est. Bit Velocity', 'Upper STD Bound', 'Lower STD Bound');


subplot(2,1,2);
% Plot error between true and estimated bit velocity
plot(t, Error_Estim_efk);
xlabel('Time [sec]');
ylabel('Error [-]','FontSize',fs)
title('Estimated Error (EKF)');