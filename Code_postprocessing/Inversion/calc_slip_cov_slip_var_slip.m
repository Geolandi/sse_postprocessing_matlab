function [slip,cov_slip,var_slip] = calc_slip_cov_slip_var_slip(L,cov_L,S,V,var_V,flag_compute_cov)

% calc_slip_cov_slip_var_slip.m computes the slip and the corresponding
% covariance matrix on the entire fault for all the epochs taken into
% account. It reconstructs the slip on the patch p at time t using the
% following equation:
% 
% slip_p(t) = sqrt(slip_strike_p(t)^2 + slip_dip_p^2)               (eq. 1)
% 
% where slip_strike_p(t) and slip_dip_p(t) are the slip on patch p at time
% t along the strike and dip directions, respectively. These quantities are
% calculated as:
% 
% slip_strike_p(t) = L_strike_p*S*V_t                              (eq. 2a)
% slip_dip_p(t)    = L_dip_p*S*V_t                                 (eq. 2b)
% 
% where L_strike and L_dip are the matrices of size PxN containing the
% result of the inversion of the spatial components U, S is a matrix of
% size NxN, and V is the matrix of size NxT. P is the total number of
% patches, N is the total number of inverted components, and T is the total
% number of epochs. Thus, L_strike_p and L_dip_p are the p-th rows of the
% L_strike and L_dip matrices, respectively, and V_t is the t-th rows of
% the V matrix.
% The function calculates also the covariance matrix of the slip at every
% single epoch using equation (S16) of the Supplementary Material of
% Gualandi et al. (2016b):
% 
%                                                                   (eq. 3)
% cov(slip_q(t),slip_r(t)) = sum_{j=1}^{3N} sum_{k=1}^{3N}
%   {
%    partiald(slip_q(t)/x_q^j(t))|_E[x_q^j(t)] * 
%    * cov(x_q^j(t),y_r^k(t)) *
%    * partiald(slip_r(t)/y_r^k(t))|_E[y_r^k(t)]
%   } =
%  = ((nabla_{x_q(t)}(slip_q(t)))|_E[x_q(t)])' * 
%    * C_{x_q(t),y_r(t)} *
%    * ((nabla_{y_r(t)}(slip_r(t)))|_E[y_r(t)])
% 
% where sum_{j=1}^{3N} and sum_{k=1}^{3N} indicate the sum from j,k=1 to
% 3N, partiald(f(z)/z) is the partial derivative of f(z) with respect to z,
% the vectors x_q(t) and y_r(t) are the set of random variables expressed
% by:
% 
% x_q(t) = [L_strike_q', L_dip_q', V_t']'                          (eq. 4a)
% y_r(t) = [L_strike_r', L_dip_r', V_t']'                          (eq. 4b)
% 
% The size of these two vectors is 3Nx1. Indeed, L_strike_p, L_dip_p, and
% V_t have size 1xN, so that L_strike_p', L_dip_p', and V_t' have size Nx1
% each. The j-th and k-th elements of these two vectors are written as
% x_q^j(t) and y_r^k(t).
% The symbol |_E[] indicates that the partial derivative is calculated in
% the point given by the expected value of the random variable indicated in
% the brackets.
% In compact matrix notation, the double sum and the elementwise partial
% derivative are written using the nabla operator for the gradient, and the
% covariance matrix C.
% Since the slip output is a matrix of size PxT, while the covariance is a
% matrix PxP for every epoch, the function outputs also the matrix with the
% variance associated to every single patch at every single epoch:
% var_slip. This matrix is useful because when we plot the time series
% relative to a given patch we use this value to estimate the standard
% deviation on the slip, i.e. its error bar, which neglects the potential
% covariances with other patches and epochs. The var_slip matrix is
% calculated using (eq. 3) for q=r=p, with p that goes from 1 to P. Using
% (eq. 3) and the assumptions of non-correlation between different
% components as well as spatial and temporal contributions, the variance of
% the slip_p(t) is given by:
% 
%                                                                   (eq. 5)
% var(slip_p(t)) = sum_{n=1}^{N}
%   [partiald(slip_p(t)/L_strike_p^n)^2 * var(L_strike_p^n) + 
%   + partiald(slip_p(t)/L_dip_p^n)^2 * var(L_dip_p^n) +
%   + 2*partiald(slip_p(t)/L_strike_p^n)*partiald(slip_p(t)/L_dip_p^n)*
%     *cov(L_strike_p^n,L_dip_p^n) +
%   + partiald(slip_p(t)/V^n(t))^2 * var(V^n(t))]
% 
% where *^n indicates the element corresponding to the n-th inverted
% component.
% -------------------------------------------------------------------------
% INPUT
% L: Matrix containing the model vectors derived inverting the spatial
%    vectors of the matrix U of the decomposition. Size: 2PxN
% cov_L: Cell variable containing the covariance matrix of the model
%        vectors. The n-th element of the cell contains the covariance
%        matrix of the n-th model vector. Cell size: Nx1. Each cell
%        contains a matrix of size: 2Px2P
% S: Diagonal matrix containing the weights associated to the inverted
%    components. Size: NxN
% V: Matrix containing the temporal information of the inverted components.
%    Size: TxN
% var_V: Matrix containing the variance associated to the matrix V. The
%        values are stored in a simple matrix because we neglect the
%        covariances through time. It is thus easier to save the variances
%        in a single matrix instead of having several diagonal matrices of
%        size TxT in a cell variable of size Nx1. Size: TxN.
% flag_compute_cov: Since the calculation of the cov_slip variable is
%                   demanding from a computation point of view, the user
%                   has the option to skip this calculation (flag must be 0
%                   in this case). If the flag is different from 0, then
%                   the function calculates the cov_slip variable. Default
%                   value: 0
% -------------------------------------------------------------------------
% OUTPUT
% slip: Matrix containing the slip value for every patch and epoch. Size:
%       PxT
% cov_slip: Cell variable containing the covariance matrix of the slip. It
%           is composed of T cells, each containing a matrix of size 2Px2P.
% var_slip: Matrix containing the variance associated to the slip value for
%           every patch and epoch. Size: PxT
% -------------------------------------------------------------------------
% References:
% Gualandi et al., Tectonophysics, 2016b. Pre- and post-seismic deformation
% related to the 2015, Mw 7.8 Gorkha Earthquake, Nepal
% 
% Adriano Gualandi - 19 Oct 2016
% California Institute of Technology
% Geological and Planetary Science Division

% If the flag_compute_cov is not given as input argument, set the default
% value to 0
if nargin<6
    flag_compute_cov = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% USEFUL VARIABLES %%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the number of patches P
P = size(L,1)/2;
% Extract the number of components inverted N and the number of epochs T
[N,T] = size(V');
% Extract the spatial slip component L relative to the strike and dip
L_strike = L(1:P,:);
L_dip    = L(P+1:2*P,:);
% Initialize the covariances and variances variables to speed up the loop
cov_L_strike_L_strike = cell(N,1);
cov_L_strike_L_dip    = cell(N,1);
cov_L_dip_L_dip       = cell(N,1);
var_L_strike = zeros(P,N);
var_L_dip    = zeros(P,N);
% For every inverted component...
for nn=1:N
    % Calculate the covariance matrices between the strike and strike
    % directions, between the strike and dip directions, and between the
    % dip and dip directions
    cov_L_strike_L_strike{nn} = cov_L{nn}(1:P,1:P);
    cov_L_strike_L_dip{nn}    = cov_L{nn}(P+1:2*P,1:P);
    cov_L_dip_L_dip{nn}       = cov_L{nn}(P+1:2*P,P+1:2*P);
    % To speed up future calculation, store the variance in a simple matrix
    % instead of being the diagonal of matrices in different cells
    var_L_strike(:,nn)  = diag(cov_L_strike_L_strike{nn});
    var_L_dip(:,nn)     = diag(cov_L_dip_L_dip{nn});
end

%%%%%%%%%%%%
%%% SLIP %%%
%%%%%%%%%%%%
% Calculate the slip along the strike and dip directions (eq. 2) and the
% final slip (eq. 1). A matrix notation is used to speed up and avoid a
% loop through the patches.
slip_strike = L_strike*S*V';
slip_dip    = L_dip*S*V';
slip = sqrt(slip_strike.^2 + slip_dip.^2);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COVARIANCE MATRIX %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% If flag_compute_cov is set to 0 then set cov_slip to NaN; otherwise
% compute the cov_slip variable
if flag_compute_cov == 0
    cov_slip     = NaN;
else
    % Initialize the variables containing the set of random variables
    % defined by (eq. 4) and the gradient of the slip with the respect to
    % these random variables used in (eq. 3).
    % See also equation (S19) of Gualandi et al. (2016b) for the definition
    % of the gradient of the slip on the patch p at time t. That equation
    % has a misprint: the last N elements should have L_strike_{p1},
    % L_dip_{p1}, ..., L_strike_{pN}, L_dip_{pN} instead of a dependence on
    % time. Furthermore, here I'm using the index n=1,...,N for the number
    % of components, while in Gualandi et al. (2016b) I used r=1,...,R. I
    % prefer now to use N to avoid confusion with the index r in the couple
    % (q,r) that goes from 1 to P.
    set_rvs = cell(P,T);
    der_slip_set_rvs = cell(P,T);
    % For every epoch...
    for tt=1:T
        % And for every patch...
        for pp=1:P
            % Define the set of random variables of (eq. 4)
            set_rvs{pp,tt} = [L_strike(pp,:), L_dip(pp,:), V(tt,:)]';
            % Find the gradient of the slip with respect to those random
            % variables (eq. (S19) of Gualandi et al., 2016b)
            der_slip_set_rvs{pp,tt} = ...
                [(S.*V(tt,:))*(slip_strike(pp,tt)/slip(pp,tt)), ...
                (S.*V(tt,:))*(slip_dip(pp,tt)/slip(pp,tt)), ...
                S.*((L_strike(pp,:)*slip_strike(pp,tt) + ...
                L_dip(pp,:)*slip_dip(pp,tt))/slip(pp,tt))]';
        end
    end
    % It remains to find the covariance matrix of the set of random
    % variables defined above. Initialize the output cov_slip variable, and
    % the variable of the covariance matrix of the set of random variables.
    cov_slip     = cell(1,T);
    cov_set_rvs  = cell(1,T);
    % Initialize also a variable to keep track of the progressing time
    elapsed_time = zeros(1,T);
    % For every epoch...
    for tt=1:T
        % Start the clock
        tic
        % Initialize the matrix cov_slip{tt} for epoch tt
        cov_slip{tt} = zeros(P,P);
        % Initialize the matrix cov_set_rvs{tt} for epoch tt
        cov_set_rvs{tt} = cell(P,P);
        % For every patch...
        for qq=1:P
            % And again for every patch...
            for rr=1:P
                % The covariance matrix that we are building here is the
                % one expressed by equation (S17) of the Supplementary
                % Material of Gualandi et al. (2016b). This is a 3Nx3N
                % matrix that must be created for every couple of patches q
                % and r, and for every epoch t. We can divide the matrix in
                % 9 blocks: top-left, top-center, top-right, middle-left,
                % middle-center, middle-right, bottom-left, bottom-center,
                % and bottom-right.
                % We initialize such 3Nx3N matrix to zeros to speed up the
                % calculations.
                % Initialize to zero the first N elements of the nn row, in
                % order to create the top-left block matrix, corresponding
                % to the covariances between L_strike_q(t) and
                % L_strike_r(t) for different inverted components nn.
                cov_set_rvs{tt}{qq,tt} = zeros(3*N,3*N);
                % TOP BLOCK MATRICES
                % For every inverted component...
                for nn=1:N
%                     cov_set_rvs{tt}{qq,rr}(nn,1:N) = zeros(1,N);
                    % Only the diagonal values of the top-left block are
                    % different from zero, and their value is given by the
                    % covariance for the component nn between L_strike_q(t)
                    % and L_strike_r(t).
                    cov_set_rvs{tt}{qq,rr}(nn,nn) = ...
                        cov_L_strike_L_strike{nn}(qq,rr);
%                     cov_set_rvs{tt}{qq,rr}(nn,N+1:2*N) = zeros(1,N);
                    % Also the top-central block matrix is diagonal, with
                    % diagonal values given by the covariance for the
                    % component nn between L_strike_q(t) and L_dip_r(t).
                    cov_set_rvs{tt}{qq,rr}(nn,N+nn) = cov_L_strike_L_dip{nn}(qq,rr);
                    % The remaining top-right matrix is zero because we
                    % assume that there is no correlation between the
                    % spatial and temporal contributions
%                     cov_set_rvs{tt}{qq,rr}(nn,2*N+1:3*N) = zeros(1,N);
                end
                % MIDDLE BLOCK MATRICES
                % For every inverted component...
                for nn=1:N
%                     cov_set_rvs{tt}{qq,rr}(N+nn,1:N) = zeros(1,N);
                    % Only the diagonal values of the central-left block
                    % are different from zero, and their value is given by
                    % the covariance for the component nn between
                    % L_strike_q(t) and L_dip_r(t).
                    cov_set_rvs{tt}{qq,rr}(N+nn,nn) = cov_L_strike_L_dip{nn}(rr,qq);
%                     cov_set_rvs{tt}{qq,rr}(N+nn,N+1:2*N) = zeros(1,N);
                    % Also the middle-central block matrix is diagonal,
                    % with diagonal values given by the covariance for the
                    % component nn between L_dip_q(t) and L_dip_r(t).
                    cov_set_rvs{tt}{qq,rr}(N+nn,N+nn) = cov_L_dip_L_dip{nn}(qq,rr);
                    % The remaining middle-right matrix is zero because we
                    % assume that there is no correlation between the
                    % spatial and temporal contributions
%                     cov_set_rvs{tt}{qq,rr}(N+nn,2*N+1:3*N) = zeros(1,N);
                end
                % BOTTOM BLOCK MATRICES
                % For every inverted component...
                for nn=1:N
                    % The only non zero block among the bottom block
                    % matrices is the bottom-right one because we are
                    % assuming that there is no correlation between the
                    % spatial and temporal contributions. The bottom-right
                    % block is diagonal because the temporal components are
                    % supposedly independent one from the other, and thus
                    % they are uncorrelated. The diagonal values are given
                    % by the variance of V(t) for the nn component.
%                     cov_set_rvs{tt}{qq,rr}(2*N+nn,1:N) = zeros(1,N);
%                     cov_set_rvs{tt}{qq,rr}(2*N+nn,N+1:2*N) = zeros(1,N);
%                     cov_set_rvs{tt}{qq,rr}(2*N+nn,2*N+1:3*N) = zeros(1,N);
                    cov_set_rvs{tt}{qq,rr}(2*N+nn,2*N+nn) = var_V(tt,nn);
                end
                % Finally, the element (q,r) of the PxP matrix cov_slip{tt}
                % is obtained using equation (S16) of the Supplementary
                % Material of Gualandi et al. (2016b), i.e. (eq. 3) in the
                % description of this function
                cov_slip{tt}(qq,rr) = ...
                    der_slip_set_rvs{qq,tt}'*cov_set_rvs{tt}{qq,rr}*...
                    der_slip_set_rvs{rr,tt};
            end
        end
        % After iterating on both qq and rr we have the full covariance
        % matrix between the slip at patch q and the slip at patch r at a
        % given time t.
%         cov_set_rvs_filename = [dirs.dir_scen,'/',dirs.dir_case,'/cov_set_rvs/cov_set_rvs_',num2str(tt),'.mat'];
%         cov_set_rvs_tt = cov_set_rvs{tt};
%         save(cov_set_rvs_filename,'cov_set_rvs_tt');
        % Keep track of the progress
        elapsed_time(tt) = toc;
        disp(['tt = ',num2str(tt),'/',num2str(T),' - Partial time: ', num2str(elapsed_time(tt)),' s,      Total time: ', num2str(sum(elapsed_time)),' s']);
    end
end

%%%%%%%%%%%%%%%%
%%% VARIANCE %%%
%%%%%%%%%%%%%%%%
% Using (eq. 5) we can iterate over the number of inverted components N,
% and using the matrix notation we can avoid the loop over the epochs T.
der_slip_L_strike = cell(N,1);
der_slip_L_dip    = cell(N,1);
der_slip_V        = cell(N,1);
var_slip_comp     = cell(N,1);
var_slip          = zeros(P,T);
for nn=1:N
    % Calculate partiald(slip_p(t)/L_strike_p^n) for all the patches and
    % epochs. This is thus a matrix PxT.
    der_slip_L_strike{nn} = ...
        (slip_strike*(S(nn,nn)).*repmat((V(:,nn))',P,1))./slip;
    % Calculate partiald(slip_p(t)/L_dip_p^n) for all the patches and
    % epochs. This is thus a matrix PxT.
    der_slip_L_dip{nn} = ...
        (slip_dip*(S(nn,nn)).*repmat((V(:,nn))',P,1))./slip;
    % Calculate partiald(slip_p(t)/V^n(t)) for all the patches and epochs.
    % This is thus a matrix PxT.
    der_slip_V{nn} = S(nn,nn)*(...
        (slip_strike.*repmat(L_strike(:,nn),1,T)) + ...
        (slip_dip.*repmat(L_dip(:,nn),1,T)))./slip;
    % Calculate the variance related to the component nn. The calculation
    % is performed simultaneously for all the patches and epochs. This is
    % thus a matrix PxT.
    var_slip_comp{nn} = ...
        (der_slip_L_strike{nn}.^2).*repmat(var_L_strike(:,nn),1,T) + ...
        (der_slip_L_dip{nn}.^2).*repmat(var_L_dip(:,nn),1,T) + ...
        (der_slip_V{nn}.^2).*repmat(var_V(:,nn)',P,1) + ...
        2*((der_slip_L_strike{nn}).*(der_slip_L_dip{nn})).*...
        repmat(diag(cov_L_strike_L_dip{nn}),1,T);
%     var_slip_comp{nn}(:,1) = var_slip_comp{nn}(:,2);
    % Calculate the final variance on the slip summing up the contributions
    % from all the inverted components.
    var_slip = var_slip + var_slip_comp{nn};
end
