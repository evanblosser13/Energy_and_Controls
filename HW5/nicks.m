% AME 5393: Renewable Energy Systems and Control
%   HW 5 : State Estimation in Oil Well Drilling
%   Nick Wandro, Sooner ID 113594933
%   Prof. Dong Zhang
%   Due May 1, 2024

% Wandro_Nick_HW5.m

clear; close all;
fs = 15;    % Font Size for plots

%% Drill String Parameters

p.J_T = 100;  % Table/top rotational inertia
p.J_B = 25;   % Bottom/bit rotational inertia
p.k = 2;      % Spring constant
p.b = 5;      % Drag coefficient

%% Problem 2 - Observability Analysis

% State space matrices
A4 = [-p.b/p.J_T, -p.k/p.J_T, 0, p.k/p.J_T;1,0,0,0;0,p.k/p.J_B, -p.b/p.J_B, -p.k/p.J_B;0, 0, 1, 0];
B4 = [1/p.J_T, 0; 0,0;0,-1/p.J_B;0,0];
C4 = [1, 0, 0, 0];

% Compute observability Matrix for 4-state system and rank
O4 = [C4;C4*A4;C4*(A4*A4);C4*(A4*A4*A4)];
disp('Rank of Observability Matrix for four-state system')
rank(O4)

% New A Matrix, for 3-state system
A = [-p.b/p.J_T, 0, -p.k/p.J_T;0,-p.b/p.J_B, p.k/p.J_T;1,-1,0];
B = [1/p.J_T, 0;0, -1/p.J_B;0,0];
C = [1,0,0];

% Observability Matrix for 3-state system and rank
O = [C;C*A;C*(A*A)];
disp('Rank of Observability Matrix for three-state system')
rank(O)

%% Problem 3 - Measurement Data
data = csvread('HW5_Data.csv');
t = data(:,1);      % t   : time vector [sec]
y_m = data(:,2);    % y_m : measured table velocity [radians/sec]
T = data(:,3);      % T   : table torque [N-m]
omega_B_true = data(:,4);    % \omega_B : true rotational speed of bit [radians/sec]

figure(1); clf;

subplot(2,1,1);
% Plot table torque
plot(t, T)
legend('Table Torque')
ylabel('[N-m]')
xlabel('Time [sec]')

subplot(2,1,2);
% Plot measured table velocity
plot(t,y_m)
legend('Table Velocity')
ylabel('[radians/sec]')
xlabel('Time [sec]')
%% Problem 4 - Luenberger Observer

% Eigenvalues of open-loop system
disp('Eigenvalues of 3-state system:')
eig(A)

% Desired poles of estimation error system
%   They should have negative real parts
%   Complex conjugate pairs
lam = [eig(A)]'*3;

% Compute observer gain (See Remark 3.1 in Notes. Use "place" command)
L = place(A',C',[6*real(lam(1));5*real(lam(1));4*real(lam(1))]);

% State-space Matrices for Luenberger Observer
A_lobs =  (A-L'*C);
B_lobs =  [1/p.J_T,L(1);0,L(2);0,0];
C_lobs =  [1, 0, 0];

sys_lobs = ss(A_lobs,B_lobs,C_lobs,[0,0]);

% Inputs to observer
u = [T, y_m];

% Initial Conditions for Luenberger Observer
x_hat0 = [0;0;0];

% Simulate Response
[y,t,x_hat] = lsim(sys_lobs,u,t,x_hat0);

% Parse out states
theta_hat = x_hat(:,3);
omega_T_hat = x_hat(:,1);
omega_B_hat = x_hat(:,2);

% Plot Results
figure(2); clf;

subplot(2,1,1);
% Plot true and estimated bit velocity
plot(t,omega_B_true,t,omega_B_hat)
legend('True','Estimated')
ylabel('Bit Speed [rad/sec]')

subplot(2,1,2);
% Plot error between true and estimated bit velocity
plot(t,(omega_B_true-omega_B_hat))
xlabel('Time [sec]')
ylabel('Error [rad/sec]')
