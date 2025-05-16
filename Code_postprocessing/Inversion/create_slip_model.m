function slip = create_slip_model(m, Cm, decomp, fault, options, ind_comps)

fprintf('Creating slip model... ')

slip.timeline = decomp.timeline;

[slip_strike,slip_dip,~] = calc_slip(...
    m,decomp.S(ind_comps,ind_comps),decomp.V(:,ind_comps));

[slip.slip,slip.cov_slip,slip.var_slip] = ...
    calc_slip_cov_slip_var_slip(m,Cm,...
    decomp.S(ind_comps,ind_comps),decomp.V(:,ind_comps),...
    decomp.var_V(:,ind_comps),options.inversion.calc_slip_cov);

%slip_dip = slip_dip - min(slip_dip,[],2);
%slip.slip = sqrt(slip_strike.^2 + slip_dip.^2);

%slip.slip(slip_dip<0) = -slip.slip(slip_dip<0);
[n_patches,n_samples] = size(slip.slip);

slip_local_vec = zeros(n_samples, n_patches, 3);
pb = ProgressBar(n_samples);
parfor t=1:n_samples
    for p=1:n_patches
        slip_local_vec(t, p, :) = fault2local(...
            [slip_strike(p,t),slip_dip(p,t),0]',...
            fault.strike(p),fault.dip(p));
    end
    pb.count;
end

slip.slip_strike = slip_strike;
slip.slip_dip = slip_dip;
slip.slip_strike_local = slip_local_vec(:,:,1)';
slip.slip_dip_local = slip_local_vec(:,:,2)';
fprintf('Done.\n')

end