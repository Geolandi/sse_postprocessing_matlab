function [U,S,V,var_U,var_V,type,llh,timeline,N,M,T] = extract_variables_from_decomp(decomp,options)

% Extract the spatial distributions U, the weights S, and the temporal
% functions V, the variances on the spatial distributions U and on the
% temporal functions V
if options.flags.flag_whitening==1
    U = decomp.DWM*decomp.U;
    var_U = (decomp.DWM.^2)*decomp.var_U;
    S = decomp.S;
    for ii=1:options.scen.N
        factU = norm(U(:,ii));
        if factU==0
            U(:,ii) = 0;
            var_U(:,ii) = Inf;
        else
            U(:,ii) = U(:,ii)/factU;
            var_U(:,ii) = var_U(:,ii)/(factU^2);
        end
        S(ii,ii) = S(ii,ii)*factU;
        clear factU;
    end
else
    U = decomp.U;
    var_U = decomp.var_U;
    S = decomp.S;
end
V = decomp.V;
var_V = decomp.var_V;

% Define the type of data used for the decomposition
type = decomp.type;
% Extract the longitude, latitude and height, and the timeline
llh = decomp.llh;
timeline = decomp.timeline;
% Check the consistency of the different quantities extracted in term of
% their sizes
[M_U,N_U] = size(U);
[T_V,N_V] = size(V);
[N_llh,~] = size(llh);
N_timeline = numel(timeline);
if N_U~=N_V
    error('The number of columns of U and V should be the same.');
else
    N = N_U;
    clear N_U N_V;
end

switch decomp.decmode
    case 'T-mode'
        if M_U~=N_llh
            error('The number of rows of U and llh should be the same.');
        else
            M = M_U;
            clear M_U N_llh;
        end
        if T_V~=N_timeline
            error(['The number of rows of V and elements of timeline ',...
                'should be the same.']);
        else
            T = T_V;
            clear T_V N_timeline;
        end
    case 'S-mode'
        if M_U~=N_timeline
            error(['The number of rows of U and elements of timeline ',...
                'should be the same.']);
        else
            M = M_U;
            clear M_U N_timeline;
        end
        if T_V~=N_llh
            error('The number of rows of V and llh should be the same.');
        else
            T = T_V;
            clear T_V N_llh;
        end
    otherwise
        error(['The decomposition mode must be one between T-mode ',...
            'and S-mode']);
end