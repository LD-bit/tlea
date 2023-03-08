function [hk_xk] = hk(xk)
% INPUT:    state vector
% OUTPUT:   observation function
x = xk(1); y = xk(2); z = xk(3);
hk_xk = [sqrt(x^2+y^2+z^2);
         acos(z/sqrt(x^2+y^2+z^2));
         atan(y/x)];
end 