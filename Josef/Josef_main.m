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

[t_cont,x_cont] = ode45(@(t, y)satelliteEOM(t, y, mu), tvec, initial_conditions, options);
y = getY(x_cont(:,1), x_cont(:,2), x_cont(:,3), x_cont(:,4), theta_0, tvec);
% plotQ1(t_cont, x_cont, y)



lx1 = 6678; % X value to linearize around
lx3 = 0; % Y value to linearize around

% df2_dx1 = -mu * (-2*lx1^2 + lx3^2) / (lx1^2 + lx3^2)^(5/2);
% df2_dx3 = mu * (3 * lx1 * lx3) / (lx1^2 + lx3^2)^(5/2);
% df4_dx1 = mu * (3 * lx1 * lx3) / (lx1^2 + lx3^2)^(5/2);
% df4_dx3 = -mu * (lx1^2 - 2*lx3^2) / (lx1^2 + lx3^2)^(5/2);

% df2_dx1 = -mu/r_0^3 + 3*mu*X_0^2/r_0^5;
% df2_dx3 = -3 * mu * X_0 * Y_0/r_0^5;
% df4_dx1 = -3 * mu * X_0 * Y_0/r_0^5;
% df4_dx3 = -mu * (1/r_0^3 - 3*Y_0^2/r_0^5);

% A_bar = [0, 1, 0, 0;...
%         -df2_dx1, 0, df2_dx3, 0;...
%         0, 0, 0, 1;...
%         df4_dx1, 0, df4_dx3, 0];

% A_bar = [0, 1, 0, 0;...
%          -mu/r_0^3, 0, 0, 0;...
%          0, 0, 0, 1;...
%          0, 0, -mu/r_0^3, 0];
temp=((X_0^2)+(Y_0^2))^(5/2);
A_barf=@(X_0,Y_0) [0 1 0 0;(2*mu*(X_0^2)-mu*(Y_0)^2)/temp 0 (3*mu*X_0*Y_0)/temp 0;0 0 0 1;(3*mu*X_0*Y_0)/temp 0 (2*mu*(Y_0^2)-mu*(X_0)^2)/temp 0];

A_bar= [0 1 0 0;(2*mu*(X_0^2)-mu*(Y_0)^2)/temp 0 (3*mu*X_0*Y_0)/temp 0;0 0 0 1;(3*mu*X_0*Y_0)/temp 0 (2*mu*(Y_0^2)-mu*(X_0)^2)/temp 0];
Fpert=eye(4)+dt.*A_bar;
% eig_A = eig(A_bar);  % Eigenvalues of the matrix A
% disp('Eigenvalues of A:');
% disp(eig_A);

B_bar = [0, 0;...
        1, 0;...
        0, 0;...
        0, 1];
A_hat = [A_bar, B_bar; zeros(2,6)];
combo = expm(A_hat * dt);
F = combo(1:4,1:4);
G = combo(1:4,5:6);

x_perturb = nan(4, length(tvec));
x_perturb(:,1) = perturb_x0;

test = nan(4, length(tvec));
test(:,1) = perturb_x0;

x_disc = nan(4, length(tvec));
x_disc(:,1) = initial_conditions;% + perturb_x0;

for i = 1:length(tvec)-1
    Fpert=eye(4)+dt.*A_barf(x_cont(i,1),x_cont(i,3));
    x_perturb(:,i+1) = Fpert * x_perturb(:, i);
    test(:,i+1) = A_bar * test(:, i);
    x_disc(:,i+1) = F * x_disc(:, i) ;%+ x_perturb(:, i);
    % x_disc(:,i+1) = F * x_disc(:, i) + x_perturb(:, i);
end

y_disc = getY(x_disc(1,:)', x_disc(2,:)', x_disc(3,:)', x_disc(4,:)', theta_0, tvec);
plotQ1(t_cont, x_cont, y)

% figure;
% subplot(4,1,1)
% plot(tvec, x_disc(1,:))
% 
% subplot(4,1,2)
% plot(tvec, x_disc(2,:))
% 
% subplot(4,1,3)
% plot(tvec, x_disc(3,:))
% 
% subplot(4,1,4)
% plot(tvec, x_disc(4,:))
% 
% 
% figure;
% subplot(4,1,1)
% plot(tvec, x_cont(:,1)+x_perturb(1,:)')
% 
% subplot(4,1,2)
% plot(tvec, x_cont(:,2)+x_perturb(2,:)')
% 
% subplot(4,1,3)
% plot(tvec, x_cont(:,3)+x_perturb(3,:)')
% 
% subplot(4,1,4)
% plot(tvec, x_cont(:,4)+x_perturb(4,:)')

figure;
subplot(4,1,1)
plot(tvec, x_perturb(1,:))

subplot(4,1,2)
plot(tvec, x_perturb(2,:))

subplot(4,1,3)
plot(tvec, x_perturb(3,:))

subplot(4,1,4)
plot(tvec, x_perturb(4,:))
% 
% figure;
% subplot(4,1,1)
% plot(tvec, test(1,:))
% 
% subplot(4,1,2)
% plot(tvec, test(2,:))
% 
% subplot(4,1,3)
% plot(tvec, test(3,:))
% 
% subplot(4,1,4)
% plot(tvec, test(4,:))







