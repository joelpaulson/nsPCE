function Jvec = map_objective_func(x, y, invCy, x0, invCx, model, y_perturb, x_perturb)

% INPUTS
% x - current parameter values (N by nx)
% y - measured outputs, i.e., data (ny by N)
% invCy - inverse of covariance of measurement noise (ny by ny) (same is used for all N)
% x0 - parameter prior mean (1 by nx)
% invCx - inverse of parameter prior covariance matrix (nx by nx)
% model - model function to evaluate (should be function x)
% y_perturb - random perturbation to data (ny by N)
% x_perturb - random perturbation to parameter (nx by 1)

% OUTPUTS
% J - objective function (scalar)

% Absolute value of x to avoid negative numbers
x = abs(x);

% Evaluate the forward model
y_pred = model(x);

% Get number
N = size(x,1);

% Transpose parameter vectors
x = x';
x0 = x0';
x_perturb = x_perturb';

% Loop over samples
Jvec = zeros(N,1);
for n = 1:N
    % Calculate the log likelihood
    J = 0;
    for i = 1:size(y,2)
        J = J + 1/2*(y(:,i) + y_perturb(:,i) - y_pred(:,i,n))'*invCy*(y(:,i) + y_perturb(:,i) - y_pred(:,i,n));
    end
    
    % Add the prior weight
    J = J + 1/2*(x(:,n) + x_perturb - x0)'*invCx*(x(:,n) + x_perturb - x0);
    
    % Store
    Jvec(n) = J;
end
end