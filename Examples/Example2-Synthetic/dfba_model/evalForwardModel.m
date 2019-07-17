function [S,T] = evalForwardModel(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that returns dfba state (S,T) for random parameter samples X
% 
% NOTE -- USER MUST MODIFY TO SPECIFY LOCAL DFBA SOLVER
% 
%         The results presented in the paper were generated using
%         DFBAlab. DFBAlab can be obtained free-of-charge for academic
%         purposes by email the authors (https://yoric.mit.edu/software/dfbalab/how-obtain-dfbalab)
%
%         Other OpenSource toolboxes exist for integrating dfba systems
%         including cobra (https://opencobra.github.io/cobratoolbox/latest/installation.html)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:   X --    matrix (N by nparam), each column should contain a list
%                   of N values of a given parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:  S --    cell array (N by 1), each element should be a matrix of
%                   dfba states evaluated over a time horizon that containts t.
%           T --    cell array (N by 1), each element should be a vector of
%                   time points corresponding to the values of S.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load initialization parameters from specified file
load('ToyProblem_dfba_parameters','model','INFO','Y0','options','tspan')

% Loop over number of samples
N = size(X,1);
INFO_nom = INFO;
model_nom = model;
Y0_nom = Y0;
S = cell(N,1);
T = cell(N,1);
print_opt = 1;
if print_opt == 1
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
    fprintf('GETTING %g SAMPLES\n', N)
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
end
for n = 1:N
    % Print statement at start of loop
    if print_opt == 1
        fprintf('integrating sample %g...', n)
        start_time = tic;
    end
    
    % Update INFO with the current sample for the parameter values
    INFO = INFO_nom;
    model = model_nom;
    Y0 = Y0_nom;
    INFO.param = X(n,:);
    Y0(1) = X(n,9);
    Y0(2) = X(n,10);
    Y0(3) = X(n,11);
    Y0(4) = X(n,12);
    model{1}.S(1,7) = -X(n,13);
    model{1}.S(2,7) = -X(n,14);
    model{1}.S(4,7) = -X(n,15);
    model{1}.S(4,4) = X(n,16);
    model{1}.S(4,5) = X(n,17);
    model{1}.S(4,6) = -X(n,18);
    
    % Setup the model
    [model,INFO] = ModelSetupM(model,Y0,INFO);
    
    % Get Lexicographic solution to LP
    if INFO.LPsolver == 0
        [INFO] = LexicographicOpt(model,INFO);
    elseif INFO.LPsolver == 1
        [INFO] = LexicographicOptG(model,INFO);
    else
        display('Solver not currently supported.');
    end
    
    % Pass to function to integrate over time for given Psamp
    [T{n},S{n},~] = intDFBA(model, tspan, Y0, options, INFO);
    
    % Print end statement
    if print_opt == 1
        end_time = toc(start_time);
        fprintf('took %g seconds\n',end_time)
    end
end
if print_opt == 1
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n')
end
end