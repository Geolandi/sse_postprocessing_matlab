function [ind_comps_common] = select_common_modes(decomp, options)
    %common_mode_perc = options["common_mode_perc"]
    common_mode_stddist = options.inversion.select_comps.common_mode_stddist;
    
    U = decomp.U;
    [n_ts, n_comp] = size(U);
    me = zeros(3, n_comp);
    st = zeros(3, n_comp);
    for j=1:3
        me(j,:) = mean(U(j:3:end,:), 1);
        st(j,:) = std(U(j:3:end,:), 1);
    end
    ind_idx_common = sum((abs(me) - common_mode_stddist*st)>0, 1);
    ind_comps_common = find(ind_idx_common==1);
end