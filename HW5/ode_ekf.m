function [xdot] = ode_efk(t,x, input_data, k1, k2, b, J_T, J_B, N, W)
    % Parse States
    theta_hat   = x(1);
    omega_T_hat = x(2);
    omega_B_hat = x(3);
    Sig = reshape(x(4:end), 3, 3);

    % Interpolate input signal data
    it = input_data(:,1);
    iy_m = input_data(:,2);
    iT = input_data(:,3);
    T= interp1(it, iT, x(4));
    y_m = interp1(it, iy_m, x(4));

    % Compute Jacobians
    F = [0, 1, -1;
        (-k1/J_T)-(3*k2/J_T)*theta_hat^2, -b/J_T, 0;
        (k1/J_B)+(3*k2/J_B)*theta_hat^2, 0, -b/J_B];
    H = [0, 1, 0];

    % Compute Kalman Gain
    L = (Sig * H' * (1/N));

    % Compute EKF system matrices
    y_hat = omega_T_hat;

    theta_hat_dot = (omega_T_hat - omega_B_hat + L(1) * (y_m - y_hat));
    omega_T_hat_dot = ((T/J_T) - (b*omega_T_hat/J_T) - (k1*theta_hat+k2*theta_hat^3)/J_T + L(2) * (y_m - y_hat));
    omega_B_hat_dot = (-(b*omega_B_hat/J_B) + (k1*theta_hat+k2*theta_hat^3)/J_B + L(3) * (y_m - y_hat));

    % Riccati Equation
    Sig_dot = ((Sig*F') + (F*Sig) + W - Sig * H' * (1/N) * H * Sig);

    % Concatenate LHS
    xdot = [theta_hat_dot; omega_T_hat_dot; omega_B_hat_dot; reshape(Sig_dot, 9, 1)];
end