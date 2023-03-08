function [Jfk_xk] = Jfk(xk,dt)
% Evaluation model jacobian matrix
% INPUT:    state vector
% OUTPUT:   transition Jacobian matrix

% VARIABLES 
Om  = 7.292115e-5;      %[rad/s] rotation rate of Earth
G   = 6.6743015e-11;    %[m^3/(kg.s^2)] gravitational constant
M   = 5.9722e24;        %[kg] Earth mass

x = xk(1); y = xk(2); z = xk(3); vx = xk(4); vy = xk(5); vz = xk(6);
% initialisation:
Jfk_xk = zeros(6,6);
% derivatives of pos function w.r.t position variables:
f1x = Om^2 + G*M*(2*x^2+y^2+z^2); 
f1y = 3*G*M*x*y;
f1z = 3*G*M*x*z;
f2x = f1y; 
f2y = Om^2 + G*M*(x^2+2*y^2+z^2);
f2z = 3*G*M*y*z;
f3x = -f1z; 
f3y = -f2z;
f3z = G*M*(x^2+y^2-z^2);
div = (x^2+y^2+z^2)^(-5/2);
Jfk_xk(1:3,1:3) = Om^2*[1 0 0;0 1 0;0 0 0] + div*[f1x f1y f1z;
                                                  f2x f2y f2z;
                                                  f3x f3y f3z];
% derivatives of pos function w.r.t velocity variables:
Jfk_xk(1:3,4:6) = 2*Om*[0 1 0;-1 0 0;0 0 0];
% derivatives of vel function w.r.t velocity variables:
Jfk_xk(4:6,4:6) = eye(3);
Jfk_xk = eye(6) + Jfk_xk*dt;
end