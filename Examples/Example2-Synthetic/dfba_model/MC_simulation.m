
%% Clear the workspace
clear
uqlab
rng(1) % fix the random seed to get repeatable results


%% Uncertainty description

% Number of MC simulations
N = 100;

% Sampling strategy for ED ['MC', 'Sobol', 'LHS', 'Halton']
%can update to other using uqlab; see their documentation for more details
sample_type = 'MC';

% Distribution of parameters
%can update to others using uqlab; see their documentation for more details
param_distribution = 'Uniform';

% Model parameter bounds
%fraction +/-
delta = 0.075;
%substrate uptake kinetics
vmaxc_lb = 1.5*(1-delta);
vmaxc_ub = 1.5*(1+delta);
Kc_lb = 0.05*(1-delta);
Kc_ub = 0.05*(1+delta);
Kie_lb = 15*(1-delta);
Kie_ub = 15*(1+delta);
vmaxn_lb = 0.25*(1-delta);
vmaxn_ub = 0.25*(1+delta);
Kn_lb = 0.5*(1-delta);
Kn_ub = 0.5*(1+delta);
vmaxo_lb = 2*(1-delta);
vmaxo_ub = 2*(1+delta);
Ko_lb = 1.2*(1-delta);
Ko_ub = 1.2*(1+delta);
vATPm_lb = 0.18*(1-delta);
vATPm_ub = 0.18*(1+delta);
%initial conditions
X0_lb = 0.01*(1-delta);
X0_ub = 0.01*(1+delta);
C0_lb = 15*(1-delta);
C0_ub = 15*(1+delta);
N0_lb = 0.3*(1-delta);
N0_ub = 0.3*(1+delta);
O0_lb = 1.0*(1-delta);
O0_ub = 1.0*(1+delta);
%coefficients for biomass creation
vxC_lb = 4*(1-delta);
vxC_ub = 4*(1+delta);
vxN_lb = 0.5*(1-delta);
vxN_ub = 0.5*(1+delta);
vxATP_lb = 1.5*(1-delta);
vxATP_ub = 1.5*(1+delta);
%coefficients for atp production
vcoxATP_lb = 1.0*(1-delta);
vcoxATP_ub = 1.0*(1+delta);
vfermATP_lb = 1.0*(1-delta);
vfermATP_ub = 1.0*(1+delta);
vlipATP_lb = 2.0*(1-delta);
vlipATP_ub = 2.0*(1+delta);
%coefficients for oxidation products
vCOXCOX_lb = 1.0*(1-delta);
vCOXCOX_ub = 1.0*(1+delta);
vfermCOX_lb = 2.0*(1-delta);
vfermCOX_ub = 2.0*(1+delta);
%stack into vector
P_lb = [vmaxc_lb, Kc_lb, Kie_lb, vmaxn_lb, Kn_lb, vmaxo_lb, Ko_lb, vATPm_lb, X0_lb, C0_lb, N0_lb, O0_lb, vxC_lb, vxN_lb, vxATP_lb, vcoxATP_lb, vfermATP_lb, vlipATP_lb, vCOXCOX_lb, vfermCOX_lb];
P_ub = [vmaxc_ub, Kc_ub, Kie_ub, vmaxn_ub, Kn_ub, vmaxo_ub, Ko_ub, vATPm_ub, X0_ub, C0_ub, N0_ub, O0_ub, vxC_ub, vxN_ub, vxATP_ub, vcoxATP_ub, vfermATP_ub, vlipATP_ub, vCOXCOX_ub, vfermCOX_ub];


%% Nominal model details

% number of models
INFO.nmodel = 1;

% stoic mat   vc vn vo vox vferm vlip vx vatp
model{1}.S = [1, 0, 0, -1, -4, -4, -4, 0 ; 
              0, 1, 0, 0, 0, 0, -0.5, 0 ;
              0, 0, 1, -1, 0, 0, 0, 0 ; 
              0, 0, 0, 1, 1, -2, -1.5, -1];

% right hand side
model{1}.b = zeros(4,1);

% the objective function 
model{1}.c = zeros(8,1);
model{1}.c(7) = 1;

% upper and lower bounds on fluxes
model{1}.lb = zeros(8,1);
model{1}.ub = inf*ones(8,1);

DB(1) = inf;
INFO.DB = DB;

% Exchange reactions
exID{1}=[1, 2, 3, 8];
INFO.exID = exID;

% Variables to decide whether to min or max
minim = 1;
maxim = -1;

% Maximize growth
C{1}(1).sense = maxim;
C{1}(1).rxns = [7];
C{1}(1).wts = [1];
% Maximize lipid
C{1}(2).sense = maxim;
C{1}(2).rxns = [6];
C{1}(2).wts = [1];
% Maximize fermentation
C{1}(3).sense = maxim;
C{1}(3).rxns = [5];
C{1}(3).wts = [1];
% Minimize carbon
C{1}(4).sense = minim;
C{1}(4).rxns = [1];
C{1}(4).wts = [1];
% Minimize nitrogen
C{1}(5).sense = minim;
C{1}(5).rxns = [2];
C{1}(5).wts = [1];
% Minimize oxygen
C{1}(6).sense = minim;
C{1}(6).rxns = [3];
C{1}(6).wts = [1];
% Minimize oxidation products
C{1}(7).sense = minim;
C{1}(7).rxns = [4];
C{1}(7).wts = [1];

% Add to INFO
INFO.C = C;

% Initial conditions
% Y1 = Biomass EColi (gDW/L)
% Y2 = Carbon (mmol/L)
% Y3 = Nitrogen (mmol/L)
% Y4 = Oxygen (mmol/L)
% Y5 = Lipids (mmol/L)
% Y6 = Ethanol (mmol/L)
% Y7 = Oxidation product (mmol/L)
% Y8 = Penalty
Y0 = [0.01, 15, 0.3, 1.0, 0.0, 0.0, 0.0, 0.0]';

% Time of simulation
tspan = [0,40];

% CPLEX Objects construction parameters
INFO.LPsolver = 0; % CPLEX = 0, Gurobi = 1.
                   % CPLEX works equally fine with both methods.
                   % Gurobi seems to work better with Method = 1, and 
                   % Mosek with Method = 0.
INFO.tol = 1E-9; % Feasibility, optimality and convergence tolerance for Cplex (tol>=1E-9). 
                 % It is recommended it is at least 2 orders of magnitude
                 % tighter than the integrator tolerance. 
                 % If problems with infeasibility messages, tighten this
                 % tolerance.
INFO.tolPh1 = INFO.tol; % Tolerance to determine if a solution to phaseI equals zero.
                   % It is recommended to be the same as INFO.tol. 
INFO.tolevt = 2*INFO.tol; % Tolerance for event detection. Has to be greater 
                   % than INFO.tol.

% You can modify the integration tolerances here.
% If some of the flows become negative after running the simulation once
% you can add the 'Nonnegative' option.
NN = 1:length(Y0);
options = odeset('AbsTol',1E-6,'RelTol',1E-6,'Nonnegative',NN,'Events',@evts);


%% Run Monte Carlo simulations

% Specify marginals
nparam = length(P_lb);
for i = 1:nparam
    Input.Marginals(i).Type = param_distribution;
    Input.Marginals(i).Parameters = [P_lb(i), P_ub(i)];
end

% Create and add the resulting input object to UQLab
myInput = uq_createInput(Input);

% Get samples from parameter distribution 
Prand = uq_getSample(myInput, N, sample_type);

% Integrate DFBA model over N samples
Y_mc = cell(N,1);
T_mc = cell(N,1);
INFO_nom = INFO;
model_nom = model;
Y0_nom = Y0;
for n = 1:N
    % Print statement at start of loop
    fprintf('integrating sample %g...', n)
    start_time = tic;
    
    % Update INFO with the current sample for the parameter values
    INFO = INFO_nom;
    model = model_nom;
    Y0 = Y0_nom;
    INFO.param = Prand(n,:);
    Y0(1) = Prand(n,9);
    Y0(2) = Prand(n,10);
    Y0(3) = Prand(n,11);
    Y0(4) = Prand(n,12);
    model{1}.S(1,7) = -Prand(n,13);
    model{1}.S(2,7) = -Prand(n,14);
    model{1}.S(4,7) = -Prand(n,15);
    model{1}.S(4,4) = Prand(n,16);
    model{1}.S(4,5) = Prand(n,17);
    model{1}.S(4,6) = -Prand(n,18);
    
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
    [Tn,Yn,INFO] = intDFBA(model, tspan, Y0, options, INFO);
    
    % Store data and print statement
    Y_mc{n} = Yn;
    T_mc{n} = Tn;
    end_time = toc(start_time);
    fprintf('took %g seconds\n',end_time)
end

% Plot results
figure
hold on
for i = 1:N
    plot(T_mc{i},Y_mc{i}(:,1),'-b',T_mc{i},Y_mc{i}(:,5),'-r','linewidth',1.5)
end
set(gcf,'color','w')
set(gca,'FontSize',16)
xlabel('time (hours)')
ylabel('concentration (g/L)')
legend('biomass', 'lipids')

figure
hold on
for i = 1:N
    plot(T_mc{i},Y_mc{i}(:,2),'-b',T_mc{i},Y_mc{i}(:,3),'-r',T_mc{i},Y_mc{i}(:,4),'-g',T_mc{i},Y_mc{i}(:,6),'-k',T_mc{i},Y_mc{i}(:,7),'-m','linewidth',1.5)
end
set(gcf,'color','w')
set(gca,'FontSize',16)
xlabel('time (hours)')
ylabel('concentration (g/L)')
legend('carbon', 'nitrogen', 'oxygen', 'ethanol', 'oxidation')

figure
hold on
for i = 1:N
    plot(T_mc{i},Y_mc{i}(:,8),'-b','linewidth',1.5)
end
set(gcf,'color','w')
set(gca,'FontSize',16)
xlabel('time (hours)')
ylabel('penalty')
