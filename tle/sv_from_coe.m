% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
function [r, v] = sv_from_coe(coe)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% This function computes the state vector (r,v) from the
% classical orbital elements (coe).
% Ref: Curtis, 2005, appendix D.9, pp. 610-611
%
% mu - gravitational parameter (kmˆ3; sˆ2)
% coe - orbital elements [h e RA incl w TA]
% where:
%       h   = angular momentum (kmˆ2/s)
%       e   = eccentricity
%       RA  = right ascension of the ascending node (rad)
%       incl = inclination of the orbit (rad)
%       w   = argument of perigee (rad)
%       TA  = true anomaly (rad)
%       R3_w - Rotation matrix about the z-axis through the angle w
%       R1_i - Rotation matrix about the x-axis through the angle i
%       R3_W - Rotation matrix about the z-axis through the angle RA
%       Q_pX - Matrix of the transformation from perifocal to
%               geocentric equatorial frame
%       rp - position vector in the perifocal frame (km)
%       vp - velocity vector in the perifocal frame (km/s)
%       r - position vector in the geocentric equatorial frame (km)
%       
%       v - velocity vector in the geocentric equatorial frame (km/s)
% 
%
% User M-functions required: none
% ------------------------------------------------------------
global mu
h   = coe(1);
e   = coe(2);
RA  = coe(3);
incl = coe(4);
w   = coe(5);
TA  = coe(6);

%...Equations 4.37 and 4.38 (rp and vp are column vectors):
rp  = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp  = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);
%...Equation 4.39:
R3_W = [cos(RA) sin(RA) 0; -sin(RA) cos(RA); 0 0 0 1];
%...Equation 4.40:
R1_i = [1 0 0; 0 cos(incl) sin(incl);0 -sin(incl) cos(incl)];
%...Equation 4.41:
R3_w = [cos(w) sin(w) 0;-sin(w) cos(w) 0;0 0 1];
%...Equation 4.44:
Q_pX = R3_W'*R1_i'*R3_w';
%...Equations 4.46 (r and v are column vectors):
r = Q_pX*rp;
v = Q_pX*vp;
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
end