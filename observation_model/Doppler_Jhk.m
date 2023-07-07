function [Jhk_xk] = Doppler_Jhk(xk,xt)
% Observation model jacobian matrix
% INPUT:    xk: state vector predicted in ECEF
%           xt: true state vector in ECEF
% OUTPUT:   observation Jacobian matrix J_{h_k}
global s_1 s_2 s_3
c = 3e8;
x = xk(1); y = xk(2); z = xk(3); vx = xk(4); vy = xk(5); vz = xk(6);
%% COEFFICIENTS

% Doppler normalised
Jh4 = (-2.1e9/c)*(vx^2+vy^2+vz^2)^(-1/2)*[vx vy vz]; 

% --> Doppler LOS (with 3 doppler measurements)
% First observer O1, expressed in ECEF: IZN-1
% Range
Jh1_1 = ((x-s_1(1))^2+(y-s_1(2))^2+(z-s_1(3))^2)^(-3/2)*[x*(x-s_1(1)), y*(y-s_1(2)), z*(z-s_1(3))]; 
% Azimuth
Jh2_1 = (sqrt((x-s_1(1))^2+(y-s_1(2))^2)*((x-s_1(1))^2+(y-s_1(2))^2+(z-s_1(3))^2))^(-1)*[(x-s_1(1))*(z-s_1(3)), (y-s_1(2))*(z-s_1(3)), -(x-s_1(1))^2-(y-s_1(2))^2+(z-s_1(3))];
% Elevation
Jh3_1 = ((y-s_1(2))^2+(x-s_1(1))^2)^(-1)*[-(y-s_1(2)), (x-s_1(1)), 0];
% Doppler LOS
Jh4_d1 = -(1/norm(xt(1:3)-s_1))*(xt(1:3)-s_1);
% Second observer O2, expressed in ECEF: note, no range at ESOC-1
% Angles
Jh2_2 = (sqrt((x-s_2(1))^2+(y-s_2(2))^2)*((x-s_2(1))^2+(y-s_2(2))^2+(z-s_2(3))^2))^(-1)*[(x-s_2(1))*(z-s_2(3)), (y-s_2(2))*(z-s_2(3)), -(x-s_2(1))^2-(y-s_2(2))^2+(z-s_2(3))];
Jh3_2 = ((y-s_2(2))^2+(x-s_2(1))^2)^(-1)*[-(y-s_2(2)), (x-s_2(1)), 0];
% Doppler LOS
Jh4_d2 = -(1/norm(xt(1:3)-s_2))*(xt(1:3)-s_2);
% Third observer O3, expressed in ECEF: note, no range at TU Graz
% Angles
Jh2_3 = (sqrt((x-s_3(1))^2+(y-s_3(2))^2)*((x-s_3(1))^2+(y-s_3(2))^2+(z-s_3(3))^2))^(-1)*[(x-s_3(1))*(z-s_3(3)), (y-s_3(2))*(z-s_3(3)), -(x-s_3(1))^2-(y-s_3(2))^2+(z-s_3(3))];
Jh3_3 = ((y-s_3(2))^2+(x-s_3(1))^2)^(-1)*[-(y-s_3(2)), (x-s_3(1)), 0];
% Doppler LOS
Jh4_d3 = -(1/norm(xt(1:3)-s_3))*(xt(1:3)-s_3);

%% FINAL JACOBIAN MATRIX
% Comment in/out the function selected below for tests:
% --> Contribution 1: Original modified function: normalised Doppler
% Jhk_xk = [Jh1_1 0 0 0;
%           Jh2_1 0 0 0;
%           Jh3_1 0 0 0;
%           0 0 0 Jh4];
% --> Contribution 2: Function with normalised Doppler + velocities 3-axis 
% Jhk_xk = [Jh1_1 0 0 0;
%           Jh2_1 0 0 0;
%           Jh3_1 0 0 0;
%           0 0 0 Jh4;
%           0 0 0 1 0 0;
%           0 0 0 0 1 0;
%           0 0 0 0 0 1];
% --> Contribution 3: Doppler with LOS
Jhk_xk = [Jh1_1 0 0 0;
          Jh2_1 0 0 0;
          Jh3_1 0 0 0;
          0 0 0 Jh4_d1';
          0 0 0 Jh4_d2';
          0 0 0 Jh4_d3'];
% --> Contribution 4: Doppler with LOS + velocity 3-axis
% Jhk_xk = [Jh1_1 0 0 0;
%           Jh2_1 0 0 0;
%           Jh3_1 0 0 0;
%           0 0 0 Jh4_d1';
%           0 0 0 1 0 0;
%           0 0 0 0 1 0;
%           0 0 0 0 0 1];
end