function l = likelihood_func(y,x,k,std_factor,std_floor,myNS,T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns the likelihood function value at a given time T(k).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   y -- vector (nq, 1), measured data at given time
%           x -- matrix (N, nparam), collection of parameter samples
%           k -- scalar, time index
%           std_factor -- vector (nq, 1), parameter of std function
%           std_floor -- vector (nq, 1), parameter of std function
%           myNS -- cell array (nq by ntime), each element is a surrogate model
%           T -- vector (ntime by 1), vector of time values that the surrogates correspond to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  l -- vector (N by 1), likelihood for all parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ny = length(y);
N = size(x,1);
mean = evalSurrogateModel(x, myNS, k, T);
std = zeros(N,ny);
for i = 1:ny
    std(:,i) = std_factor(i)*mean(:,i)+std_floor(i);
end
l = prod(normpdf(repmat(y,[N,1]),mean,std),2);
end