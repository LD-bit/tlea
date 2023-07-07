function [Rot] = rotationENU_ECEF(theta,phi)
% Create the rotation matrix from ENU to ECEF with measurements of 
% spherical ENU coordinates azimuth and elevation
% INPUTS:   theta, azimuth [rad]
%           phi, elevation [rad]
Rot = [-sin(theta) -cos(theta)*sin(phi) cos(phi)*cos(theta);
       cos(theta) -sin(phi)*sin(theta) cos(phi)*sin(theta);
       0 cos(phi) sin(phi)];
end