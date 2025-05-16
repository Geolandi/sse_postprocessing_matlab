function G = create_greens_function_tri_NikkhooandWalter2015(stations, fault, nu, origin, flag_parallel)
    if ~exist('nu','var')
        nu = 0.25;
    end
    if ~exist('origin','var')
        origin = [mean(mean(fault.lon)), mean(mean(fault.lat))];
    end
    
    km2m = 1e3;
    
    n_patches = size(fault.lon,1);
    
    v1.lon = zeros(n_patches,1);
    v1.lat = zeros(n_patches,1);
    v1.height = zeros(n_patches,1);
    v2.lon = zeros(n_patches,1);
    v2.lat = zeros(n_patches,1);
    v2.height = zeros(n_patches,1);
    v3.lon = zeros(n_patches,1);
    v3.lat = zeros(n_patches,1);
    v3.height = zeros(n_patches,1);
    ind_reorganize = zeros(n_patches,3);
    for j=1:n_patches
        [~, ind_lon] = sort(fault.lon(j,:),'ascend');
        ind1 = ind_lon(1);
        ind23 = [1,2,3];
        ind23(ind1) = [];
        [~, ind_lat] = sort(fault.lat(j,ind23),'ascend');
        ind_reorganize(j,1) = ind1;
        ind_reorganize(j,2:3) = ind23(ind_lat);
        v1.lon(j,:)    = fault.lon(j,ind_reorganize(j,1));
        v1.lat(j,:)    = fault.lat(j,ind_reorganize(j,1));
        v1.height(j,:) = fault.height(j,ind_reorganize(j,1));
        v2.lon(j,:)    = fault.lon(j,ind_reorganize(j,2));
        v2.lat(j,:)    = fault.lat(j,ind_reorganize(j,2));
        v2.height(j,:) = fault.height(j,ind_reorganize(j,2));
        v3.lon(j,:)    = fault.lon(j,ind_reorganize(j,3));
        v3.lat(j,:)    = fault.lat(j,ind_reorganize(j,3));
        v3.height(j,:) = fault.height(j,ind_reorganize(j,3));
    end
    
    % v1.lon = fault.lon(:,1);
    % v1.lat = fault.lat(:,1);
    % v1.height = fault.height(:,1);
    % 
    % v2.lon = fault.lon(:,2);
    % v2.lat = fault.lat(:,2);
    % v2.height = fault.height(:,2);
    % 
    % v3.lon = fault.lon(:,3);
    % v3.lat = fault.lat(:,3);
    % v3.height = fault.height(:,3);


    v1 = llh2localxyz_geodetic(v1, origin);
    v2 = llh2localxyz_geodetic(v2, origin);
    v3 = llh2localxyz_geodetic(v3, origin);
    
    stations.origin = origin;
    stations = llh2localxyz_geodetic(stations, origin);

    n_patches = length(v1.lon);
    n_stations = length(stations.lon);
    
    tic
    fprintf('Creating Greens'' function:\n')
    G = zeros(n_stations*3, n_patches*2);
    if flag_parallel > 1
        % Construct a ParforProgressbar object
        pb = ProgressBar(n_patches);
        % For every patch...
        Gstrike = zeros(n_stations*3, n_patches);
        Gdip = zeros(n_stations*3, n_patches);
        parfor j = 1:n_patches
            Gj = zeros(n_stations*3, 2);
    
            x = [v1.xE(j), v2.xE(j), v3.xE(j)]*km2m;
            y = [v1.yN(j), v2.yN(j), v3.yN(j)]*km2m;
            z = [v1.zV(j), v2.zV(j), v3.zV(j)]*km2m;
            for i = 1:n_stations
                sx = stations.xE(i)*km2m;
                sy = stations.yN(i)*km2m;
                sz = 0; %stations.zV(i)*km2m;
                
                [ue_ss,un_ss,uv_ss] = TDdispHS(sx,sy,sz,...
                    [x(1), y(1), z(1)],...
                    [x(2), y(2), z(2)],...
                    [x(3), y(3), z(3)],1,0,0,nu);
                [ue_ds,un_ds,uv_ds] = TDdispHS(sx,sy,sz,...
                    [x(1), y(1), z(1)],...
                    [x(2), y(2), z(2)],...
                    [x(3), y(3), z(3)],0,1,0,nu);
                
                Gj(i*3-2, 1) = ue_ss;
                Gj(i*3-1, 1) = un_ss;
                Gj(i*3, 1)   = uv_ss;
    
                Gj(i*3-2, 2) = ue_ds;
                Gj(i*3-1, 2) = un_ds;
                Gj(i*3, 2)   = uv_ds;
            end
            Gstrike(:,j) = Gj(:,1);
            Gdip(:,j) = Gj(:,2);
            pb.count;
        end
        G(:,1:n_patches) = Gstrike;
        G(:,n_patches+1:2*n_patches) = Gdip;
    else
        for j = 1:n_patches
            if mod(j,500) == 0
              fprintf('Patch %d/%d: %.2f s\n', j, n_patches, toc)
            end
            x = [v1.xE(j), v2.xE(j), v3.xE(j)]*km2m;
            y = [v1.yN(j), v2.yN(j), v3.yN(j)]*km2m;
            z = [v1.zV(j), v2.zV(j), v3.zV(j)]*km2m;
            for i = 1:n_stations
                sx = stations.xE(i)*km2m;
                sy = stations.yN(i)*km2m;
                sz = 0; %stations.zV(i)*km2m;
                
                [ue_ss,un_ss,uv_ss] = TDdispHS(sx,sy,sz,...
                    [x(1), y(1), z(1)],...
                    [x(2), y(2), z(2)],...
                    [x(3), y(3), z(3)],1,0,0,nu);
                [ue_ds,un_ds,uv_ds] = TDdispHS(sx,sy,sz,...
                    [x(1), y(1), z(1)],...
                    [x(2), y(2), z(2)],...
                    [x(3), y(3), z(3)],0,1,0,nu);
                
                G(i*3-2, j) = ue_ss;
                G(i*3-1, j) = un_ss;
                G(i*3, j)   = uv_ss;

                G(i*3-2, n_patches+j) = ue_ds;
                G(i*3-1, n_patches+j) = un_ds;
                G(i*3, n_patches+j)   = uv_ds;
            end
        end
    end
    fprintf('Done in %.2f s\n', toc)
    
end