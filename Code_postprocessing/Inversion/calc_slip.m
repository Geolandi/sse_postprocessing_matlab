function [slip_strike,slip_dip,slip] = calc_slip(L,S,V)

% calc_slip.m calculates the slip distribution at P patches for T epochs
% using the decomposition of N components.
% -------------------------------------------------------------------------
% INPUT
% L: Matrix containing the model vectors derived inverting the spatial
%    vectors of the matrix U of the decomposition. Size: 2PxN
% S: Diagonal matrix containing the weights associated to the inverted
%    components. Size: NxN
% V: Matrix containing the temporal information of the inverted components.
%    Size: TxN
% -------------------------------------------------------------------------
% OUTPUT
% slip_strike: Matrix containing the slip value for every patch and epoch
%              in the strike direction. Size: PxT
% slip_dip: Matrix containing the slip value for every patch and epoch in
%           the dip direction. Size: PxT
% slip: Matrix containing the slip value for every patch and epoch. Size:
%       PxT
% -------------------------------------------------------------------------
% 
% Adriano Gualandi - 20 Oct 2016
% California Institute of Technology
% Geological and Planetary Sciences Division

% Extract the number of patches P
P = size(L,1)/2;
% Extract the spatial slip component L relative to the strike and dip
L_strike = L(1:P,:);
L_dip    = L(P+1:2*P,:);
% Calculate the slip along the strike and dip directions (eq. 2) and the
% final slip (eq. 1). A matrix notation is used to speed up and avoid a
% loop through the patches.
slip_strike = L_strike*S*V';
slip_dip    = L_dip*S*V';
slip = sqrt(slip_strike.^2 + slip_dip.^2);
