function [X, Y] = getTrackingStationPos(omega_e, theta_0, R_e, time_vec)
% Author: Josef Michelsen
% Date: 11/26/2024

num_stations = length(theta_0);
num_timesteps = length(time_vec);

X = nan(num_stations, num_timesteps);
Y = nan(num_stations, num_timesteps);

for i = 1:num_stations
    for j = 1:num_timesteps
        t = time_vec(j);
        X(i, j) = R_e * cos(omega_e * t + theta_0(i));
        Y(i, j) = R_e * sin(omega_e * t + theta_0(i));
    end
end

end

