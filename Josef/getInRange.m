function in_range = getInRange(theta, phi)
% Author: Josef Michelsen
% Date: 12/6/2024

upper_bound = wrapToPi(pi/2 + theta);
lower_bound = wrapToPi(-pi/2 + theta);

in_range = false;
if upper_bound > lower_bound
    in_range = phi > lower_bound && phi < upper_bound;
elseif upper_bound < lower_bound
    in_range = phi < upper_bound || phi > lower_bound;
end

end

