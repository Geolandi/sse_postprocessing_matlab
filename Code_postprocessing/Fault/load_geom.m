function [fault, G] = load_geom(decomp, options, dirs)

% load fault matrix
fault_matrix = load([dirs.dir_raid,dirs.dir_data,...
        options.fault.fault_file]);
% each row corresponds to a patch; for each patch there are 3 vertices,
% and thus 3 values of lon, 3 values of lat, and 3 values of height
fault.lon = [fault_matrix(:,1), fault_matrix(:,4), fault_matrix(:,7)];
fault.lat = [fault_matrix(:,2), fault_matrix(:,5), fault_matrix(:,8)];
fault.height = [fault_matrix(:,3), fault_matrix(:,6), fault_matrix(:,9)];
% each row of llh is the lon, lat and height of the center of each patch
fault.llh = [mean(fault.lon,2), mean(fault.lat,2), mean(fault.height,2)];
clear fault_matrix

% total number of patches
n_patches = length(fault.lon(:,1));

%% CALCULATE GREENS' FUNCTIONS
% locate the origin at the center of the mesh
origin = [mean(mean(fault.lon)), mean(mean(fault.lat))];
fault.origin = origin;
% transform geografic (degrees) coordinates into local (km) coordinates
fault.xE = zeros(n_patches,3);
fault.yN = zeros(n_patches,3);
fault.zV = zeros(n_patches,3);
for i=1:3
    fault_i.lon = fault.lon(:,i);
    fault_i.lat = fault.lat(:,i);
    fault_i.height = fault.height(:,i);
    fault_i = llh2localxyz_geodetic(fault_i, origin);
    fault.xE(:,i) = fault_i.xE;
    fault.yN(:,i) = fault_i.yN;
    fault.zV(:,i) = fault_i.zV;
end

fault.xyz = [mean(fault.xE,2), mean(fault.yN,2), mean(fault.zV,2)];

switch decomp.name{1}(end)
    case 'u'
        data_type = '1d';
    case 'e'
        switch decomp.name{3}(end)
            case 'e'
                data_type = '2d';
            case 'u'
                data_type = '3d';
            otherwise

        end
    otherwise

end

switch data_type
    case '1d'
        llh_stations = decomp.llh(1:end,:);
    case '2d'
        llh_stations = decomp.llh(1:2:end,:);
    case '3d'
        llh_stations = decomp.llh(1:3:end,:);
    otherwise

end

G = create_greens_function(llh_stations, fault, options, origin);

%% Find strike and rake of patches
fault.strike = zeros(n_patches,1);
fault.dip    = zeros(n_patches,1);
fault.area   = zeros(n_patches,1);
for i=1:n_patches
    % A: x, y, and z coordinates of the first point
    % B: x, y, and z coordinates of the second point
    % C: x, y, and z coordinates of the third point
    A = [fault.xE(i,1), fault.yN(i,1), fault.zV(i,1)];
    B = [fault.xE(i,2), fault.yN(i,2), fault.zV(i,2)];
    C = [fault.xE(i,3), fault.yN(i,3), fault.zV(i,3)];
    [a,b,c,d] = plane_3points(A,B,C);
    fault.strike(i) = 90 - atan2d(-a, b);
    % atan2d returns values between -180 and 180, the strike is thus
    % between -90 and 270
    % if -a/c<0, atan2d must be between 0 and 180, so the strike must be
    % between -90 and 90
    if -a/c < 0
        if atan2d(-a, b) < 0
            fault.strike(i) = fault.strike(i) + 180;
        end
    else
        if atan2d(-a, b) > 0
            fault.strike(i) = fault.strike(i) - 180;
        end
    end
    % let us keep the strike between 0 and 360
    fault.strike = mod(fault.strike, 360);
    fault.dip(i) = 90 - atan2d(c, sqrt(a^2 + b^2));
    if fault.dip(i)>90
        fault.dip(i) = 180 - fault.dip(i);
    end
    AB = norm(B-A);
    AC = norm(C-A);
    theta = acosd(dot(B-A, C-A) / (AB * AC));
    fault.area(i) = 0.5*AB*AC*sind(theta);
end

% fault2.xE = fault.xyz(:,1);
% fault2.yN = fault.xyz(:,2);
% fault2.zV = fault.xyz(:,3);
% fault2.strike = fault.strike;
% fault2.dip = fault.dip;
% [n,P_fault,strike_fault,dip_fault] = norm_plane(fault2);

end