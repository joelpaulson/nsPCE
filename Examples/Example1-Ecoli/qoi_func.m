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
Y = zeros(N,3);

% Loop over number of samples
for i = 1:N
    % Get unique time points for sample i
    [T_interp,ia,~] = unique(T{i}');
    
    % Interpolate for glucose value at t
    Yg_interp = S{i}(ia,3)';
    Y(i,1) = interp1(T_interp,Yg_interp,t);

    % Interpolate for xylose value at t
    Yx_interp = S{i}(ia,4)';
    Y(i,2) = interp1(T_interp,Yx_interp,t);

    % Interpolate for biomass value at t
    Yb_interp = S{i}(ia,2)';
    Y(i,3) = interp1(T_interp,Yb_interp,t);
end
end