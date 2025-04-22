function [xdot] = ode_kf(t,x,A,B,C,input_data,W,N)

% Parse States
x_hat = x(1:3);
Sig = reshape(x(4:end),[3 3]);   % Need to reshape Sigma from vector to matrix

% Parse and interpolate input signal data
%   This is a coding trick for time-varying input data.
it = input_data(:,1);

iy_m = input_data(:,2);
iT = input_data(:,3);

T = interp1(it,iT,t);
y_m = interp1(it,iy_m,t);

% Compute Kalman Gain (Look at Chapter 3, Section 4)
L = Sig * C' * (1/N);

% Kalman Filter equations
x_hat_dot = A * x_hat + B * T + L * (y_m - (C * x_hat));

% Riccati Equation for Sigma
Sig_dot = Sig * A' + A * Sig + W - Sig * C' * (1/N) * C * Sig;

% Concatenate LHS
xdot = [x_hat_dot; reshape(Sig_dot, [9 1])];