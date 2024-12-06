function y = getYLinear(X, theta_0, time_vec)
% Author: Josef Michelsen
% Date: 12/6/2024
if size(X, 1) ~= 4
    X = X';
end

[X_s, X_s_d, Y_s, Y_s_d] = getTrackingStationStates(theta_0, time_vec);

rho = @(x, xs, y, ys) sqrt((x - xs)^2 + (y - ys)^2);

c = @(x,xs,xd,xds,y,ys,yd,yds) [(x-xs)/rho(x,xs,y,ys), 0, (y-ys)/rho(x,xs,y,ys), 0; ...
    ((xd-xds)*rho(x,xs,y,ys)-(x-xs)*((x-xs)*(xd-xds)+(y-ys)*(yd-yds)))/rho(x,xs,y,ys)^3, (x-xs)/rho(x,xs,y,ys), ((yd-yds)*rho(x,xs,y,ys)-(y-ys)*((x-xs)*(xd-xds)+(y-ys)*(yd-yds)))/rho(x,xs,y,ys)^3, (y-ys)/rho(x,xs,y,ys);...
    (y-ys)/rho(x,xs,y,ys)^2, 0, (x-xs)/rho(x,xs,y,ys)^2,0];

% c = @(x,xs,xd,xds,y,ys,yd,yds) [(x-xs)/rho(x,xs,y,ys), 0, (y-ys)/rho(x,xs,y,ys), 0; ...
%     ((y-ys)*((y-ys)*(xd-xds)-(x-xs)*(yd-yds)))/rho(x,xs,y,ys)^3, (x-xs)/rho(x,xs,y,ys), ((x-xs)*((x-xs)*(yd-yds)-(y-ys)*(xd-xds)))/rho(x,xs,y,ys)^3, (y-ys)/rho(x,xs,y,ys);...
%     (y-ys)/rho(x,xs,y,ys), 0, (x-xs)/rho(x,xs,y,ys),0];

y = cell(length(theta_0), length(X));

for i = 1:length(theta_0)
    for j = 1:length(X)
        C = c(X(1,j), X_s(i,j), X(2,j), X_s_d(i,j), X(3,j), Y_s(i,j), X(4,j), Y_s_d(i,j));
        y_temp = C * X(:,j);
        in_range = getInRange(y_temp(3), y_temp(2));

        if in_range
            y{i,j} = y_temp;
        else
            y{i,j} = nan(3,1);
        end
    end
end

end

