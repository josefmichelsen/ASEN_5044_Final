% Author: Josef Michelsen
% Date: 11/26/2024

clear; close all; clc;

load("../orbitdeterm_finalproj_KFdata.mat");

% Constants
mu = 398600; % Standard gravitational parameter [km^3/s^2]
R_e = 6378; % Radius of the Earth [km]
omega_e = 2 * pi / 86400; % Turning rate of Earth [rad/s]
dt = 10;

% Initial Conditions
X_0 = 6678; % [km]
r_0 = X_0; % [km]
Y_0 = 0; % [km]
X_d_0 = 0; % [km/s]
Y_d_0 = sqrt(mu/r_0); % [km/s]
initial_conditions = [X_0; X_d_0; Y_0; Y_d_0];

% Initial ground station positions
stations = 1:12;
theta_0 = (stations - 1) .* pi/6;
% X_s_0 = R_e .* cos(theta_0);
% Y_s_0 = R_e .* sin(theta_0);
% 
% [X_t, Y_t] = getTrackingStationPos(omega_e, theta_0, R_e, tvec);

perturb_x0 = [0; 0.075; 0; -0.021];
initial_conditions_cont = initial_conditions + perturb_x0;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

[t_cont, x_cont] = ode45(@(t, y)satelliteEOM(t, y, mu), tvec, initial_conditions_cont, options);
y = getY(x_cont(:,1), x_cont(:,2), x_cont(:,3), x_cont(:,4), theta_0, tvec);
% plotQ1(t_cont, x_cont, y)

[t_cont_noPerturb, x_cont_noPerturb] = ode45(@(t, y)satelliteEOM(t, y, mu), tvec, initial_conditions, options);



lx1 = 6678; % X value to linearize around
lx3 = 0; % Y value to linearize around

denom = @(X_0,Y_0) ((X_0^2)+(Y_0^2))^(5/2);

A_bar_func = @(X_0,Y_0) [0, 1, 0, 0; ...
                    (2*mu*(X_0^2)-mu*(Y_0)^2)/denom(X_0,Y_0), 0, (3*mu*X_0*Y_0)/denom(X_0,Y_0), 0; ...
                    0, 0, 0, 1; ...
                    (3*mu*X_0*Y_0)/denom(X_0,Y_0), 0, (2*mu*(Y_0^2)-mu*(X_0)^2)/denom(X_0,Y_0), 0];

A_bar = A_bar_func(X_0, Y_0);

Fpert = eye(4) + dt .* A_bar;

B_bar = [0, 0;...
        1, 0;...
        0, 0;...
        0, 1];

x_perturb = nan(4, length(tvec));
x_perturb(:,1) = perturb_x0;

% x_disc = nan(4, length(tvec));
% x_disc(:,1) = initial_conditions;% + perturb_x0;

for i = 1:length(tvec)-1
    Fpert = eye(4) + dt .* A_bar_func(x_cont_noPerturb(i,1), x_cont_noPerturb(i,3));
    x_perturb(:,i+1) = Fpert * x_perturb(:, i);

    % x_disc(:,i+1) = Fpert * x_disc(:, i);
    % x_disc(:,i+1) = F * x_disc(:, i) + x_perturb(:, i);
end

x_combo = x_cont_noPerturb' + x_perturb;

y_linear_disc = getYLinear(x_perturb, x_cont_noPerturb ,theta_0, tvec, y);
% y_linear_disc = getYLinear(x_combo, x_cont_noPerturb ,theta_0, tvec);
% y_disc = getY(x_combo(1,:)', x_combo(2,:)', x_combo(3,:)', x_combo(4,:)', theta_0, tvec);
% plotQ1(t_cont, x_combo, y_disc)
plotQ1(t_cont, x_combo, y_linear_disc)

% plotStates(tvec, x_perturb)







