function slip_rate = create_slip_rate_model(slip, fault, options)

fprintf('Creating slip rate model... ')
[n_patches,n_samples] = size(slip.slip);
slip_rate_strike = zeros(n_patches, n_samples)*NaN;
slip_rate_dip = zeros(n_patches, n_samples)*NaN;
slip_strike = slip.slip_strike;
slip_dip = slip.slip_dip;
timeline = slip.timeline;

slip_rate_method = options.slip_rate.method;
switch slip_rate_method
    case 'tvdiff'
        tvdiff_iter = options.slip_rate.tvdiff.iter;
        tvdiff_alph = options.slip_rate.tvdiff.alph;
        tvdiff_scale = options.slip_rate.tvdiff.scale;
        tvdiff_ep = options.slip_rate.tvdiff.ep;
        tvdiff_dx_string = options.slip_rate.tvdiff.dx;
        switch tvdiff_dx_string
            case 'dt'
                tvdiff_dx = mean(diff(timeline));
            case '1d'
                tvdiff_dx = 1/365.25;
            otherwise
                error('Unrecognized `options.slip_rate.tvdiff.dx`.')
        end
        tvdiff_plotflag = options.slip_rate.tvdiff.plotflag;
        tvdiff_diagflag = options.slip_rate.tvdiff.diagflag;
        rate_window = NaN;
    case 'slopemovwin'
        rate_window = options.slip_rate.slopemovwin.win_length;
        tvdiff_iter = NaN;
        tvdiff_alph = NaN;
        tvdiff_scale = NaN;
        tvdiff_ep = NaN;
        tvdiff_dx = NaN;
        tvdiff_plotflag = NaN;
        tvdiff_diagflag = NaN;
    case 'slopemovwincen'
        rate_window = options.slip_rate.slopemovwincen.win_length;
        tvdiff_iter = NaN;
        tvdiff_alph = NaN;
        tvdiff_scale = NaN;
        tvdiff_ep = NaN;
        tvdiff_dx = NaN;
        tvdiff_plotflag = NaN;
        tvdiff_diagflag = NaN;
    otherwise

end


% Construct a ParforProgressbar object
pb = ProgressBar(n_patches);
parfor j=1:n_patches
    %fprintf('%d\n', j)
    slip_strike_j = slip_strike(j,:);
    slip_dip_j = slip_dip(j,:);
    switch slip_rate_method
        case 'tvdiff'
            v_strike = TVRegDiff(slip_strike_j,...
                tvdiff_iter,tvdiff_alph,[],tvdiff_scale,tvdiff_ep,...
                tvdiff_dx,tvdiff_plotflag,tvdiff_diagflag)';
            v_dip = TVRegDiff(slip_dip_j,...
                tvdiff_iter,tvdiff_alph,[],tvdiff_scale,tvdiff_ep,...
                tvdiff_dx,tvdiff_plotflag,tvdiff_diagflag)';
            slip_rate_strike(j,:) = v_strike(end-n_samples+1:end);
            slip_rate_dip(j,:) = v_dip(end-n_samples+1:end);
        case 'slopemovwin'
            slip_rate_strike_j = zeros(1, n_samples) * NaN;
            slip_rate_dip_j = zeros(1, n_samples) * NaN;
            for t=rate_window:n_samples
                p_strike = polyfit(...
                    timeline(t-rate_window+1:t)-timeline(t-rate_window+1), ...
                    slip_strike_j(t-rate_window+1:t), 1);
                p_dip = polyfit(...
                    timeline(t-rate_window+1:t)-timeline(t-rate_window+1), ...
                    slip_dip_j(t-rate_window+1:t), 1);
                slip_rate_strike_j(t) = p_strike(1);
                slip_rate_dip_j(t) = p_dip(1);
            end
            slip_rate_strike(j,:) = slip_rate_strike_j;
            slip_rate_dip(j,:) = slip_rate_dip_j;
        case 'slopemovwincen'
            slip_rate_strike_j = zeros(1, n_samples) * NaN;
            slip_rate_dip_j = zeros(1, n_samples) * NaN;
            half_rate_window = floor(0.5*rate_window);
            for t=half_rate_window+1:n_samples-half_rate_window
                t1 = t-half_rate_window;
                t2 = t+half_rate_window;
                p_strike = polyfit(...
                    timeline(t1:t2)-timeline(t1), ...
                    slip_strike_j(t1:t2), 1);
                p_dip = polyfit(...
                    timeline(t1:t2)-timeline(t1), ...
                    slip_dip_j(t1:t2), 1);
                slip_rate_strike_j(t) = p_strike(1);
                slip_rate_dip_j(t) = p_dip(1);
            end
            slip_rate_strike(j,:) = slip_rate_strike_j;
            slip_rate_dip(j,:) = slip_rate_dip_j;
        otherwise

    end
    pb.count;
end

rake_requested = options.slip_rate.rake;
rake_lim1 = rake_requested - 90;
rake_lim2 = rake_requested + 90;
rake_obs = atan2d(slip_rate_dip,slip_rate_strike);
ind_flip = rake_obs < rake_lim1 | rake_obs > rake_lim2;
slip_rate.slip_rate = sqrt(slip_rate_strike.^2 + slip_rate_dip.^2);
slip_rate.slip_rate(ind_flip) = -slip_rate.slip_rate(ind_flip);
% slip_rate.slip_rate(slip_rate_dip<0) = ...
%     -slip_rate.slip_rate(slip_rate_dip<0);
slip_rate.slip_rate_strike = slip_rate_strike;
slip_rate.slip_rate_dip = slip_rate_dip;
fprintf('Done.\n')

fprintf('Converting slip rate vectors for plots... ')
slip_rate_local_vec = zeros(n_samples, n_patches, 3);
pb = ProgressBar(n_samples);
parfor t=1:n_samples
    for p=1:n_patches
        slip_rate_local_vec(t, p, :) = fault2local(...
            [slip_rate_strike(p,t),slip_rate_dip(p,t),0]',...
            fault.strike(p),fault.dip(p));
    end
    pb.count;
end
slip_rate.slip_rate_strike_local = slip_rate_local_vec(:,:,1)';
slip_rate.slip_rate_dip_local = slip_rate_local_vec(:,:,2)';
fprintf('Done.\n')
end