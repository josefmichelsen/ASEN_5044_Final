function [X, Xd, Y, Yd] = getTrackingStationStates(theta_0, time_vec)
% Author: Josef Michelsen
% Date: 11/26/2024

R_e = 6378; % Radius of the Earth [km]
omega_e = 2 * pi / 86400; % Turning rate of Earth [rad/s]

num_stations = length(theta_0);
num_timesteps = length(time_vec);

X = nan(num_stations, num_timesteps);
Xd = nan(num_stations, num_timesteps);
Y = nan(num_stations, num_timesteps);
Yd = nan(num_stations, num_timesteps);

for i = 1:num_stations
    for j = 1:num_timesteps
        t = time_vec(j);
        theta = theta_0(i);

        X(i, j) = R_e * cos(omega_e * t + theta);
        Xd(i, j) = -omega_e * R_e * sin(omega_e * t + theta);
        Y(i, j) = R_e * sin(omega_e * t + theta);
        Yd(i,j) = omega_e * R_e * cos(omega_e * t + theta);
    end
end

end

