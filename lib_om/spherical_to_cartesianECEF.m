function [x] = spherical_to_cartesianECEF(lat,lon,alt)
% Conversion of spherical lat,lon,alt coordinates into Cartesian (x,y,z)
% in ECEF reference frame.
% INPUTS:   latitude [degrees]
%           longitude [degrees]
%           altitude above sea level [m]

% VARIABLES
R_e_equator = 6378.137e3;           %[m] Earth's radius equator
R_e_poles = 6356.7523e3;            %[m] Earth's radius poles
R_e =  (R_e_poles+R_e_equator)/2;   %[m] Earth's radius average

x = (R_e+alt)*[cosd(lat)*cosd(lon);
               cosd(lat)*sind(lon);
               sind(lat)];
end