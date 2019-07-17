function [Y] = qoi_func(S,T,t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the quantities of interest Y at a specific time
% point t given samples of the dfba states (S,T). The state information
% should be evaluated over a time horizon, and this data will be
% interpolated to determine the quantity of interest at any chosen t. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   S --    cell array (N by 1), each element should be a matrix of
%                   dfba states evaluated over a time horizon that containts t.
%           T --    cell array (N by 1), each element should be a vector of
%                   time points corresponding to the values of S.
%           t --    scalar, the time point of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  Y --    matrix (N by nq), each column corresponds to a specific
%                   quantity of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract number of samples
N = length(S);

% Initialize the output matrix
Y = zeros(N,7);

% Loop over number of samples
for i = 1:N
    % Get unique time points for sample i
    [T_interp,ia,~] = unique(T{i}');

    % Interpolate for biomass value at t
    Y(i,1) = interp1(T_interp,S{i}(ia,1)',t);    
    
    % Interpolate for carbon value at t
    Y(i,2) = interp1(T_interp,S{i}(ia,2)',t);

    % Interpolate for nitrogen value at t
    Y(i,3) = interp1(T_interp,S{i}(ia,3)',t);

    % Interpolate for oxygen value at t
    Y(i,4) = interp1(T_interp,S{i}(ia,4)',t);

    % Interpolate for lipid value at t
    Y(i,5) = interp1(T_interp,S{i}(ia,5)',t);

    % Interpolate for ethanol value at t
    Y(i,6) = interp1(T_interp,S{i}(ia,6)',t);

    % Interpolate for oxidation value at t
    Y(i,7) = interp1(T_interp,S{i}(ia,7)',t);    
end
end