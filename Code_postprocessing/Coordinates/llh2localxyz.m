function [x,y,z] = llh2localxyz(lat,lon,height,origin)

% llh2localxyz.m converts the geographic coordinates of a list of points
% into local coordinates.
% -------------------------------------------------------------------------
% INPUT
% lat: latitude (degrees)
% lon: longitude (degrees)
% height: altitude (km)
% origin: [longitude, latitude] (degrees) of the origin of the new local
%         coordinate reference system
% -------------------------------------------------------------------------
% OUTPUT
% x: x in the local reference system (km), i.e. along the East direction
% y: y in the local reference system (km), i.e. along the North direction
% z: z in the local reference system (km), i.e. along the Vertical
%    direction
% -------------------------------------------------------------------------
%
% Adriano Gualandi - 19 Aug 2016
% California Institute of Technology
% Geological and Planetary Science Division

xy = (llh2localxy([lat, lon, height]',fliplr(origin)));
x = xy(:,1);
y = xy(:,2);
z = height;