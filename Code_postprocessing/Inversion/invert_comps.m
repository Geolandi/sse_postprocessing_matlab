function [m, Cm] = invert_comps(ICA, ind_comps, fault, G, ind_sigma0_comps, options)

%%
n_ICs2invert = length(ind_comps);

n_patches = length(fault.lon(:,1));

%m0 = zeros(2*n_patches,1);


%lambda_strike = options.inversion.lambda_strike;
%lambda_dip = options.inversion.lambda_dip;
% d_patches_triu = triu(d_patches); min(d_patches_triu(d_patches_triu>0))
%lambda0 = options.inversion.lambda0;
% lambda_strike = 2*sqrt(mean(fault.area));
% lambda_dip = 2*sqrt(mean(fault.area));
% lambda0 = 2*sqrt(mean(fault.area));
% etm0s = options.inversion.sigma(ind_sigma);
%which_smoothing = options.inversion.which_smoothing;

% n_smoothings = length(etm0s);

dobs = ICA.U(:,ind_comps);
Cdobs = cell(1,n_ICs2invert);
for i=1:n_ICs2invert
    Cdobs{i} = diag(ICA.var_U(:,ind_comps(i)));
end

%m = cell(1,n_smoothings);
%Cm = cell(1,n_smoothings);

% compute distances between patches in the local coordinates
d_patches = pdist2(fault.xyz(:,1:3),fault.xyz(:,1:3));

if options.flags.flag_parallel > 1
    % Construct a ParforProgressbar object
    pb = ProgressBar(options.scen.n_comps_max);
    
    m = zeros(2*n_patches,n_ICs2invert);
    Cm = cell(1, n_ICs2invert);
    parfor i=1:n_ICs2invert
        if options.inversion.verbose == true
            fprintf('   IC #%d\n', ind_comps(i))
        end
        [m0, Cm0] = create_priors(n_patches, options, ind_sigma0_comps(i), ...
            d_patches);
        Cm{i} = zeros(2*n_patches,2*n_patches);
        [m(:,i),Cm{i}(:,:)] = invert_comp(G,...
            dobs(:,i), ...
            Cdobs{i},...
            m0,...
            Cm0);
    end
    pb.count;
    if options.inversion.verbose == true
        fprintf('Done\n')
        fprintf('\n')
    end
    
else
    m = zeros(2*n_patches,n_ICs2invert);
    Cm = cell(1, n_ICs2invert);
    for i=1:n_ICs2invert
        if options.inversion.verbose == true
            fprintf('   IC #%d\n', ind_comps(i))
        end
        [m0, Cm0] = create_priors(n_patches, options, ind_sigma0_comps(i), ...
            d_patches);
        Cm{i} = zeros(2*n_patches,2*n_patches);
        [m(:,i),Cm{i}(:,:)] = invert_comp(G,...
            dobs(:,i), ...
            Cdobs{i},...
            m0,...
            Cm0);
    end
    %save([dirs.dir_scen,dirs.dir_case,'/matfiles/inversion_sigma0_',...
    %    num2str(etm0),'.mat'],'m','Cm','-v7.3');
    if options.inversion.verbose == true
        fprintf('Done\n')
        fprintf('\n')
    end
end

end