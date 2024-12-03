function plotQ1(t, x, y)
% Author: Josef Michelsen
% Date: 11/30/2024

figure;
subplot(4,1,1)
plot(t,x(:, 1))
xlabel("Time [sec]")
ylabel("X [km]")

subplot(4,1,2)
plot(t,x(:, 2))
xlabel("Time [sec]")
ylabel("Xdot [km]")

subplot(4,1,3)
plot(t,x(:, 3))
xlabel("Time [sec]")
ylabel("Y [km]")


subplot(4,1,4)
plot(t,x(:, 4))
xlabel("Time [sec]")
ylabel("Ydot [km]")

sgtitle("States vs. Time, Full Nonlinear Dynamics Simulation")
fontsize(16,'points')

station = cell(size(y, 1), 1);
for i = 1:size(y, 1)
temp_y = reshape([y{i,:}], 3, []);
    temp_rho = temp_y(1,:);
    temp_rho_d = temp_y(2,:);
    temp_phi = temp_y(3,:);
    station{i} = [temp_rho; temp_rho_d; temp_phi]; 
end

figure;
subplot(4,1,1)
hold on
for i = 1:size(y, 1)
    scatter(t, station{i}(1,:))
end
hold off
ylabel("Rho^i [km]")

subplot(4,1,2)
hold on
for i = 1:size(y, 1)
    scatter(t, station{i}(2,:))
end
hold off
ylabel("Rhodot^i [km/s]")

subplot(4,1,3)
hold on
for i = 1:size(y, 1)
    scatter(t, station{i}(3,:))
end
hold off
ylabel("$\Phi^i$ [rad/s]", 'Interpreter','latex')

subplot(4,1,4)
hold on
for i = 1:size(y, 1)
    visible_vec = ~isnan(station{i}(3,:)) * i;
    visible_vec(visible_vec == 0) = nan;
    scatter(t, visible_vec)
end
hold off
ylabel("Visible Station ID")
xlabel("Time [sec]")
sgtitle("Full Nonlinear Model Data Simulation")
fontsize(16,'points')

end

