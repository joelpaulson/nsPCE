function [T, Y, INFO] = intDFBA(model, tspan, Y0, options, INFO)
% Integrates DFBA model, need to update with whatever solver available

tint = 0;
TF = [];
YF = [];
while tint<tspan(2)
    % Look at MATLAB documentation if you want to change solver.
    % ode15s is more or less accurate for stiff problems.
    [T,Y] = ode15s(@DRHS,tspan,Y0,options,INFO);
    TF = [TF;T];
    YF = [YF;Y];
    tint = T(end);
    tspan = [tint,tspan(2)];
    Y0 = Y(end,:);
    if tint == tspan(2)
        break;
    end
    
    %Determine model with basis change
    value = evts(tint,Y0,INFO);
    [jjj,j] = min(value);
    ct = 0;
    k = 0;
    while j>ct
        k = k + 1;
        ct = ct + size(model{k}.A,1);
    end
    INFO.flagbasis = k;
%     fprintf('Basis change at time %d. \n',tint);
    
    % Update b vector
    [INFO] = bupdate(tint,Y0,INFO);
    
    % Perform lexicographic optimization
    if INFO.LPsolver == 0
        [INFO] = LexicographicOpt(model,INFO);
    elseif INFO.LPsolver == 1
        [INFO] = LexicographicOptG(model,INFO);
    else
        display('Solver not currently supported.');
    end
end
T = TF;
Y = YF;

end