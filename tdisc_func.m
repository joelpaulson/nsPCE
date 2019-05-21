function [Td] = tdisc_func(S,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the discontinuity time td given samples of the 
% dfba states (S,T). The state information should be evaluated over a time
% horizon. The columns should be ordered sequentially such that the first
% column corresponds to the first discontinuity, etc. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   S --    cell array (N by 1), each element should be a matrix of
%                   dfba states evaluated over a time horizon.
%           T --    cell array (N by 1), each element should be a vector of
%                   time points corresponding to the values of S.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  Td --   matrix (N by ndisc), each column containts a unique
%                   discontinuity time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract number of samples
N = length(S);

% Initialize the output
Td = zeros(N,2);

% Loop over number of samples
for i = 1:length(S)
    % Find time that glucose is zero
    index = find(S{i}(:,3)<1e-10,1);
    if ~isempty(index)
        Td(i,1) = T{i}(index);        
    else
        Td(i,1) = inf;
    end
    
    % Find time that xylose is zero
    index = find(S{i}(:,4)<1e-10,1);
    if ~isempty(index)
        Td(i,2) = T{i}(index);        
    else
        Td(i,2) = inf;
    end
end
end