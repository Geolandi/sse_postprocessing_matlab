function [geodetic_output]  = llh2localxyz_geodetic(geodetic_input,origin)

% llh2localxyz_geodetic.m add to the geodetic structure the xE, yN, and zV
% fields containing the catalog coordinate in the local reference system.
% -------------------------------------------------------------------------
% INPUT
% geodetic_input: structure that must have the following fields:
%   lat: latitude of the geodetic data (degrees)
%   lon: longitude of the geodetic data (degrees)
%   height: height of the geodetic data (km)
% origin: [longitude, latitude] (degrees) of the origin of the new local
%         coordinate reference system
% -------------------------------------------------------------------------
% OUTPUT
% geodetic_output: structure equal to geodetic_input but with the
% additional fields:
%   xE: East location of the geodetic data in the local reference system
%       (km)
%   yN: North location of the geodetic data in the local reference system
%       (km)
%   zV: Vertical location of the geodetic data in the local reference
%       system (km)
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 26 Aug 2016
% California Institute of Technology
% Geological and Planetary Science Division

geodetic_output = geodetic_input;
[geodetic_output.xE,geodetic_output.yN,geodetic_output.zV] = ...
    llh2localxyz(geodetic_input.lat,geodetic_input.lon,...
    geodetic_input.height,origin);
