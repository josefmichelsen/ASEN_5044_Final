function xdot = satelliteEOM(time, input, mu)
% Author: Josef Michelsen
% Date: 11/29/2024

X = input(1);
X_d = input(2);
Y = input(3);
Y_d = input(4);

r = sqrt(X^2 + Y^2);

X_dd = -mu * X / r^3;
Y_dd = -mu * Y / r^3;

xdot = [X_d; X_dd; Y_d; Y_dd];

end

