function [m0, Cm0] = create_priors(n_patches, options, ind_sigma, d_patches)


lambda_strike = options.inversion.lambda_strike;
lambda_dip = options.inversion.lambda_dip;
lambda0 = options.inversion.lambda0;
etm0s = options.inversion.sigma;
which_smoothing = options.inversion.which_smoothing;

m0 = options.inversion.m0*ones(2*n_patches,1);

etm0 = etm0s(ind_sigma);

if options.inversion.verbose == true
    fprintf('sigma0 = %f\n', etm0);
end

Em0 = ones(1,n_patches) * etm0;

% normalize standard deviation on model depending on the smoothing
% (see Valette, 05/11/2009)
Em_strike = Em0.*lambda0./lambda_strike;
Em_dip    = Em0.*lambda0./lambda_dip;

switch which_smoothing
    case 'gaussian' 
        % use bsxfun to multiply each line of Cm by Em
        Cm0_strike = bsxfun(@times, ...
            exp(-0.5*d_patches.^2/lambda_strike^2),Em.^2);
        Cm0_dip = bsxfun(@times, ...
            exp(-0.5*d_patches.^2/lambda_dip^2),Em.^2);
    case 'exponential'
        Cm0_strike = bsxfun(@times, ...
            exp(-d_patches/lambda_strike),Em_strike.^2);
        Cm0_dip = bsxfun(@times, ...
            exp(-d_patches/lambda_strike),Em_dip.^2);

    otherwise
        error([which_smoothing 'is not a valid option for ' ...
            'smoothing choice. Choices are '...
            '`expotential` or `gaussian`']);
end

Cm0 = [Cm0_strike, zeros(n_patches,n_patches);...
       zeros(n_patches,n_patches), Cm0_dip];

end