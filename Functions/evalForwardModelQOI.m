function y = evalForwardModelQOI(x, tlist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns all quantities of interest y for given parameter x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   x --    matrix (1 by nparam), each column should contain a list
%                   of N values of a given parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  y --    matrix (nq by ntime), each element is one quantitity of
%                   interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S,T] = evalForwardModel(x);

y = [];
for i = 1:length(tlist)
    y = [y ; qoi_func(S,T,tlist(i))];
end
y = y';

end