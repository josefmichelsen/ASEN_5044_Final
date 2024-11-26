% Author: Josef Michelsen
% Date: 11/26/2024

clear; close all; clc;

load("../orbitdeterm_finalproj_KFdata.mat");

% Constants
mu = 398600; % Standard gravitational parameter [km^3/s^2]
R_e = 6378; % Radius of the Earth [km]
omega_e = 2 * pi / 86400; % Turning rate of Earth [rad/s]

% Initial Conditions
X_0 = 6678; % [km]
r_0 = X_0; % [km]
Y_0 = 0; % [km]
X_d_0 = 0; % [km/s]
Y_d_0 = r_0 * sqrt(mu/r_0^3); % [km/s]

% Initial ground station positions
theta_0 = (1:12 - 1) .* pi/6;
X_s_0 = R_e .* cos(theta_0);
Y_s_0 = R_e .* sin(theta_0);

[X_t, Y_t] = getTrackingStationPos(omega_e, theta_0, R_e, tvec);

