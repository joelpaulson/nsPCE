function Y = evalSurrogateModel(X, myNS, k, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the evaluates a pce surrogate at a particular time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   X    -- matrix (N by nparam), each column corresponds to a parameter
%           myNS -- cell array (nq by ntime), each element is a surrogate model
%           k    -- scalar, the time index of interest \in [1, ntime]
%           T    -- vector (ntime by 1), vector of time values that the surrogates correspond to 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  Y    -- matrix (N by nq), each column corresponds to a specific
%                   quantity of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(X,1);
nq = size(myNS,1);
Y = zeros(N,nq);
for iq = 1:nq
    if strcmp(myNS{iq,k}.fit_type, 'ns')
        Td_val = uq_evalModel(myNS{iq,k}.td.myPCE, X);
        indexRegion{1} = find(Td_val > T(k));
        indexRegion{2} = find(Td_val <= T(k));
    else
        indexRegion{1} = 1:N;
    end
    for ik = 1:myNS{iq,k}.Ne
        Y(indexRegion{ik},iq) = uq_evalModel(myNS{iq,k}.mesh{ik}.myPCE, X(indexRegion{ik},:));
    end
    index_neg = find(Y(:,iq) < 0);
    Y(index_neg,iq) = zeros(length(index_neg),1);    
end
end