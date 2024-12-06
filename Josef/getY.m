function y = getY(X, X_d, Y, Y_d, theta_0, time_vec)
% Author: Josef Michelsen
% Date: 11/30/2024

% R_e = 6378; % Radius of the Earth [km]
% omega_e = 2 * pi / 86400; % Turning rate of Earth [rad/s]

[X_s, X_s_d, Y_s, Y_s_d] = getTrackingStationStates(theta_0, time_vec);

% figure;
% hold on
% for i = 1:size(X_s, 1)
%     scatter3(X_s(i,:), Y_s(i,:), repmat(i, size(X_s(i,:))))
% end
% hold off

rho = nan(length(theta_0), length(X));
rho_d = nan(length(theta_0), length(X));
phi = nan(length(theta_0), length(X));
theta_i = nan(length(theta_0), length(X));
y = cell(length(theta_0), length(X));

for i = 1:length(theta_0)
    for j = 1:length(X)
        rho(i,j) = sqrt((X(j) - X_s(i,j))^2 + (Y(j) - Y_s(i,j))^2);

        % X_s_d = -omega_e * R_e * sin(omega_e * time_vec(j) + theta_0(i));
        % Y_s_d = omega_e * R_e * cos(omega_e * time_vec(j) + theta_0(i));

        rho_d(i,j) = ((X(j) - X_s(i,j)) * (X_d(j) - X_s_d(i,j)) + (Y(j) - Y_s(i,j)) * (Y_d(j) - Y_s_d(i,j))) / rho(i,j);
        
        phi(i,j) = atan2(Y(j) - Y_s(i,j), X(j) - X_s(i,j));
        
        theta_i(i,j) = atan2(Y_s(i,j), X_s(i,j));

        % upper_bound = wrapToPi(pi/2 + theta_i(i,j));
        % lower_bound = wrapToPi(-pi/2 + theta_i(i,j));
        % 
        % in_range = false;
        % if upper_bound > lower_bound
        %     in_range = phi(i,j) > lower_bound && phi(i,j) < upper_bound;
        % elseif upper_bound < lower_bound
        %     in_range = phi(i,j) < upper_bound || phi(i,j) > lower_bound;
        % end
        in_range = getInRange(theta_i(i,j), phi(i,j));

        if in_range
            y{i,j} = [rho(i,j); rho_d(i,j); phi(i,j)];
        else
            y{i,j} = nan(3,1);
        end
    end
end

end

