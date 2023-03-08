function [fk_xk] = fk(xk)
% Evaluation model function
% INPUT:    state vector
% OUTPUT:   transition function

% VARIABLES 
Om  = 7.292115e-5;      %[rad/s] rotation rate of Earth
G   = 6.6743015e-11;    %[N.m^2/(kg.s^2)] gravitational constant
M   = 5.9722e24;        %[kg] Earth mass

x = xk(1); y = xk(2); z = xk(3); vx = xk(4); vy = xk(5); vz = xk(6);
fk_xk = [Om^2*x + 2*Om*vy - (G*M*x)/((x^2 + y^2 + z^2)^(3/2));
         Om^2*y - 2*Om*vx - (G*M*y)/((x^2 + y^2 + z^2)^(3/2));
         - (G*M*z)/((x^2 + y^2 + z^2)^(3/2));
         vx; 
         vy;
         vz];
end