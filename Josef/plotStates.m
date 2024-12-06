function plotStates(t, x)
% Author: Josef Michelsen
% Date: 12/02/2024

if size(x, 2) ~= 4
    x = x';
end

figure;
subplot(4,1,1)
plot(t, x(:, 1), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("X [km]")

subplot(4,1,2)
plot(t, x(:, 2), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Xdot [km]")

subplot(4,1,3)
plot(t, x(:, 3), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Y [km]")


subplot(4,1,4)
plot(t, x(:, 4), 'LineWidth', 2)
xlabel("Time [sec]")
ylabel("Ydot [km]")

sgtitle("States vs. Time")
fontsize(16,'points')


end

