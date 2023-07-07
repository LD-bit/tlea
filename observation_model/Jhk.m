function [Jhk_xk] = Jhk(xk)
% Observation model jacobian matrix
% INPUT:    state vector
% OUTPUT:   observation Jacobian matrix
x = xk(1); y = xk(2); z = xk(3);
% initialisation:
Jh1 = (x^2+y^2+z^2)^(-3/2)*[x y z];
Jh2 = (sqrt(x^2+y^2)*(x^2+y^2+z^2))^(-1)*[x*z y*z -x^2-y^2+z^2];
Jh3 = (y^2+x^2)^(-1)*[-y x 0];
Jhk_xk = [Jh1 0 0 0;
          Jh2 0 0 0;
          Jh3 0 0 0];
end