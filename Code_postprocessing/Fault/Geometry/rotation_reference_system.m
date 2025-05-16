function R = rotation_reference_system(axis,angle)

% rotation_reference_system.m gives the linear operator that allows to
% perform a rotation of the reference system around a given axis of a given
% angle.
% -------------------------------------------------------------------------
% INPUT
% axis: string that can assume the value 'x', 'y', or 'z'
% angle: rotation angle in degrees, increasing counterclockwise in a right
% hand system
% -------------------------------------------------------------------------
% OUTPUT
% R: linear operator that can be applied to a 3-d vector to change the
% reference system in the rotated one
% -------------------------------------------------------------------------
% 
% Adriano Gualandi
% California Institute of Technology
% Geological and Planetary Science Division

switch axis
    case 'x'
        R = [ 1, 0, 0;...
             0, cosd(angle), sind(angle);...
             0, -sind(angle), cosd(angle)];
    case 'y'
        R = [cosd(angle), 0, -sind(angle);...
             0, 1, 0;...
             sind(angle), 0, cosd(angle)];
    case 'z'
        R = [ cosd(angle), sind(angle), 0;...
             -sind(angle), cosd(angle), 0;...
              0, 0, 1];
    otherwise
        error('Please, select a string between ''x'', ''y'', and ''z'' for the axis variable.');
end