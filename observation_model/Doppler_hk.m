function [hk_xk] = Doppler_hk(xk,xt)
% INPUT:    xk: state vector predicted in ECEF
%           xt: true state vector in ECEF
% OUTPUT:   observation function output [R, el, az, vr_1, vr_2, vr_3]
global s_1 s_2 s_3
c = 3e8;
x = xk(1); y = xk(2); z = xk(3); vx = xk(4); vy = xk(5); vz = xk(6);
%% COEFFICIENTS
% First observer O1, expressed in ECEF: IZN-1
R1          = sqrt((x-s_1(1))^2+(y-s_1(2))^2+(z-s_1(3))^2); % O1: range
theta1      = acos((z-s_1(3))/sqrt((x-s_1(1))^2+(y-s_1(2))^2+(z-s_1(3))^2));% O1: azimtuh
phi1        = atan((y-s_1(2))/(x-s_1(1))); % O1: elevation
% O1: Doppler
XusVu1      = [xt(1) - s_1(1), xt(2) - s_1(2), xt(3) - s_1(3)]*[xt(4) xt(5) xt(6)]';
fd1_norm    = -(1/(norm(xt(1:3)-s_1)))*XusVu1;

% Second observer O2, expressed in ECEF: note, no range at ESOC-1
theta2      = acos((z-s_2(3))/sqrt((x-s_2(1))^2+(y-s_2(2))^2+(z-s_2(3))^2)); % O2: azimtuh
phi2        = atan((y-s_2(2))/(x-s_2(1))); % O2: elevation
% O2: Doppler
XusVu2      = [xt(1) - s_2(1), xt(2) - s_2(2), xt(3) - s_2(3)]*[xt(4) xt(5) xt(6)]';
fd2_norm    = -(1/(norm(xt(1:3)-s_2)))*XusVu2;

% Third observer O2, expressed in ECEF: note, no range at TU Graz
theta3      = acos((z-s_3(3))/sqrt((x-s_3(1))^2+(y-s_3(2))^2+(z-s_3(3))^2)); % O3: azimtuh
phi3        = atan((y-s_3(2))/(x-s_3(1))); % O3: elevation
% O3: Doppler
XusVu3      = [xt(1) - s_3(1), xt(2) - s_3(2), xt(3) - s_3(3)]*[xt(4) xt(5) xt(6)]';
fd3_norm    = -(1/(norm(xt(1:3)-s_3)))*XusVu3;

%% FINAL OBSERVATION FUNCTION
% Comment in/out the function selected below for tests:
% --> Contribution 1: Original modified function: normalised Doppler
% hk_xk = [R1;...
%          theta1;...
%          phi1;...
%          -(2.1e9/c)*(vx^2+vy^2+vz^2)^(1/2)];
% --> Contribution 2: Function with normalised Doppler + velocities 3-axis 
% hk_xk = [R1;...
%          theta1;...
%          phi1;...
%          -(1/c)*(vx^2+vy^2+vz^2)^(1/2);...
%          vx;...
%          vy;...
%          vz];
% --> Contribution 3: Doppler with LOS
hk_xk = [R1;theta1;phi1;fd1_norm;fd2_norm;fd3_norm];
% --> Contribution 4: Doppler with LOS + velocity 3-axis
% hk_xk = [R1;theta1;phi1;fd1_norm;vx;vy;vz];
end 