function y = evalSurrogateModelFull(x, myNS, time_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the evaluates a pce surrogate for all times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   X         -- matrix (N by nparam), each column corresponds to a parameter
%           myNS      -- cell array (nq by ntime), each element is a surrogate model
%           time_list -- vector (ntime by 1), vector of time values that the surrogates correspond to 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  Y    -- matrix (nq by ntime by N), samples of all quantites of
%                   interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ny = size(myNS,1);
N = size(myNS,2);
Nsamp = size(x,1);

y = zeros(ny,N,Nsamp);
for i = 1:ny
    for j = 1:N
        if strcmp(myNS{i,j}.fit_type, 'ns')
            Td_val = uq_evalModel(myNS{i,j}.td.myPCE, x);
            indexRegion{1} = find(Td_val > time_list(j));
            indexRegion{2} = find(Td_val <= time_list(j));
        else
            indexRegion{1} = 1:Nsamp;
        end
        for ik = 1:myNS{i,j}.Ne
            y(i,j,indexRegion{ik}) = uq_evalModel(myNS{i,j}.mesh{ik}.myPCE, x(indexRegion{ik},:));
        end
        index_neg = find(y(i,j,:) < 0);
        y(i,j,index_neg) = zeros(length(index_neg),1);
    end
end

end