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
load('EColi_dfba_parameters','model','INFO','Y0','options','tspan')

% Loop over number of samples
N = size(X,1);
INFO_nom = INFO;
model_nom = model;
S = cell(N,1);
T = cell(N,1);
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('GETTING %g SAMPLES\n', N)
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
for n = 1:N
    % Print statement at start of loop
    fprintf('integrating sample %g...', n)
    start_time = tic;
    
    % Update INFO with the current sample for the parameter values
    INFO = INFO_nom;
    model = model_nom;
    INFO.param = X(n,:);
    
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
    end_time = toc(start_time);
    fprintf('took %g seconds\n',end_time)
end
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n')
end