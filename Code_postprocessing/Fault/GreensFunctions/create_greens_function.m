function G = create_greens_function(llh_stations, fault, options, origin)

stations.lon = llh_stations(:,1);
stations.lat = llh_stations(:,2);
stations.height = llh_stations(:,3);
% compute triangular Greens' functions
switch options.fault.method
    case 'Meade2007'
        G = create_greens_function_tri_Meade2007(...
            stations,fault,options.fault.nu,origin);
    case 'NikkhooandWalter2015'
        G = create_greens_function_tri_NikkhooandWalter2015(...
            stations,fault,options.fault.nu,origin,...
            options.flags.flag_parallel);
    case 'pointPCAIM'
        G = create_greens_function_pointPCAIM(...
            stations,fault,options.fault.nu,origin);
    case 'rect'
        G = create_greens_function_rect(...
            stations,fault,options.fault.nu,origin);
    otherwise

end

end