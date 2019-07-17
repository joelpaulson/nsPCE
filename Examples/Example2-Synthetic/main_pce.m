
%% Initialize

% Start uqlab -- NOTE: MUST DOWNLOAD UQLAB (free for academic users https://www.uqlab.com/)
uqlab

% Clear the workspace
clear

% Fix the random seed to get repeatable results
rng(120)



%% User inputs

% Quantity of interest function
%function of dfba states
qoi = @(S,T,t)(qoi_func(S,T,t));
%time points of interest
tlist = [10, 20, 30, 40];
%add dfba integrator to path -- NOTE: MUST BE MODIFIED BY USER (see evalForwardModel.m for more information)
addpath('./dfba_model')

% Choice of pce implementation ['ns' or 'global']
fit_type = 'ns';

% If 'ns', must specify discontinuity function
tdisc = @(S,T)(tdisc_func(S,T));

% Uncertain parameters
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
%probability distribution
param_distribution = 'Uniform';

% Algorithm parameters
%tolerance -- modified leave-one-out (LOO) cross-validated error
eLOO_tol = 1e-3;
%initial number of experimental design (ED) samples
N_init = 500;
%number of ED samples to sequentially add if tolerance is not met
N_add = 100;
%maximum number of ED samples
N_max = 2500;
%factor multiplying tolerance to decide if element is constant or not
std_factor = 1.0;
%pce parameters -- see uqlab documentation for more info
metaopts.Type = 'Metamodel';
metaopts.MetaType = 'PCE';
metaopts.Method = 'LARS' ;                  % select sparse-favoring least squares
metaopts.Degree = 1:30;                     % degree of polynomial. If array is given then automatically select one with lowest LOO error 
metaopts.TruncOptions.qNorm = 0.6;          % sparse truncation scheme
metaopts.LARS.LarsEarlyStop = true;         % specify early stop
%sampling strategy for ED ['MC', 'Sobol', 'LHS', 'Halton'] -- see uqlab for more info
sample_type = 'MC';

% Save the workspace? Empty vector = No
save_mat_file = [];

% Validation file? Empty vector = No -- NOTE: MUST DEFINE "Y_val" and "X_val" to match qoi function
validation_mat_file = [];



%% Generate the initial experimental design (ED) data set

% Specify marginal distribution
nparam = length(P_lb);
for i = 1:nparam
    Input.Marginals(i).Type = param_distribution;
    Input.Marginals(i).Parameters = [P_lb(i), P_ub(i)];
end

% Create and add the resulting input object to uqlab
myInput = uq_createInput(Input);

% Sample the parameter space
N_ED = N_init;
X_ED = uq_getSample(myInput, N_ED, sample_type);

% Integrate the full dfba model at ED points
[S_ED,T_ED] = evalForwardModel(X_ED);

% Evaluate qoi for all time points
ntime = length(tlist);
nq = size(qoi(S_ED(1),T_ED(1),tlist(1)),2);
Y_ED = zeros(N_ED,nq,ntime);
for it = 1:ntime
    Y_ED(:,:,it) = qoi(S_ED,T_ED,tlist(it));
end



%% Fit pce for every discontinuity time

if strcmp(fit_type, 'ns')
    % Calculate ED for Td
    Td_ED = tdisc(S_ED,T_ED);
    ndisc = size(Td_ED,2);
    
    % For all discontinuities, fit global pce
    tdisc_pce = cell(ndisc,1);
    for i = 1:ndisc
        % Initialize the variables
        tdisc_pce{i}.myInput = myInput;
        tdisc_pce{i}.eLOO = inf;
        
        % Run until error tolerance is met
        while tdisc_pce{i}.eLOO > eLOO_tol
            % Add samples to ED after first iteration
            if tdisc_pce{i}.eLOO ~= inf
                % Sample the parameters space
                X_ED_new = uq_getSample(tdisc_pce{i}.myInput, N_add, sample_type);                
                % Integrate the full dfba model at new samples
                [S_ED_new,T_ED_new] = evalForwardModel(X_ED_new);                
                % Evaluate the qoi for all time points and all new samples
                Y_ED_new = zeros(N_add,nq,ntime);
                for it2 = 1:ntime
                    Y_ED_new(:,:,it2) = qoi(S_ED_new,T_ED_new,tlist(it2));
                end                
                % Evaluate the discontinuity time for all new samples
                Td_ED_new = tdisc(S_ED_new,T_ED_new);                
                % Add samples to global ED
                N_ED = N_ED + N_add;
                X_ED = [X_ED ; X_ED_new];
                Y_ED = [Y_ED ; Y_ED_new];
                Td_ED = [Td_ED ; Td_ED_new];
            end            
            % Create and add to uqlab the PCE model
            metaopts.Input = myInput;
            metaopts.ExpDesign.X = X_ED;
            metaopts.ExpDesign.Y = Td_ED(:,i);
            tdisc_pce{i}.myInput = myInput;
            tdisc_pce{i}.myPCE = uq_createModel(metaopts);            
            % Extract the modified LOO error
            tdisc_pce{i}.eLOO = tdisc_pce{i}.myPCE.Error.ModifiedLOO;            
            % Break out of loop if max number of samples exceeded
            if N_ED >= N_max
                break
            end
        end
        
        % Get samples for initialization purposes when fitting pce in each element
        N_pred = 1e4;
        tdisc_pce{i}.X = uq_getSample(tdisc_pce{i}.myInput, N_pred, sample_type);
        tdisc_pce{i}.Y = uq_evalModel(tdisc_pce{i}.myPCE, tdisc_pce{i}.X);
    end
    
elseif strcmp(fit_type,'global')
    % Do nothing
    
else
    % Print error if fit_type is not set properly
    error('The fit_type must be either "ns" or "global"')
end



%% Get pce representation for every quantity of interest

% Initialize cell storing all models
myNS = cell(nq,ntime);
for iq = 1:nq
    for it = 1:ntime
        myNS{iq,it}.fit_type = fit_type;
    end
end

% Print statement
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('STARTING PCE CONSTRUCTION FOR QOI\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n')

% Loop over time
for it = 1:ntime
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If ns, check if discontinuity has occured
    index_disc = zeros(ndisc,1);
    idisc = 1;
    if strcmp(myNS{1,it}.fit_type, 'ns')
        for idisc = ndisc:-1:1
            [Ftd,xtd]=ecdf(tdisc_pce{idisc}.Y);
            index = find(Ftd > 0.05,1);
            if tlist(it) < xtd(index)
                index = [];
            else
                index = 1.0;
            end
            if ~isempty(index)
                index_disc(idisc) = 1;
                break
            end
        end
    end

    % If a discontinuity has not happened, switch to global pce
    if index_disc(idisc) ~= 1
        for iq = 1:nq
            myNS{iq,it}.fit_type = 'global';
            myNS{iq,it}.Ne = 1;
        end    
    % Otherwise add the dis
    else
        for iq = 1:nq
            myNS{iq,it}.idisc = idisc;
            myNS{iq,it}.td = tdisc_pce{idisc};
            myNS{iq,it}.Ne = 2;
        end        
    end
    
    % If ns, then initialize two mesh
    if strcmp(myNS{1,it}.fit_type, 'ns')
        % Classify point in ED for which discontinuity has not happened (class=1) and has happened (class=2)
        iClass{1} = find(Td_ED(:,idisc) > tlist(it));
        iClass{2} = find(Td_ED(:,idisc) <= tlist(it));
        
        % Estimate probability of each region
        probClass{1} = length(find(tdisc_pce{idisc}.Y > tlist(it)))/N_pred;
        probClass{2} = 1 - probClass{1};
        
        % Define two element mesh
        for iq = 1:nq
            for k = 1:myNS{iq,it}.Ne
                [aClass, bClass] = getBoundsClass(myNS{iq,it}.td, tlist(it), k, P_lb, P_ub);
                myNS{iq,it}.mesh{k}.a = aClass;
                myNS{iq,it}.mesh{k}.b = bClass;
                myNS{iq,it}.mesh{k}.N_ED = length(iClass{k});
                myNS{iq,it}.mesh{k}.X_ED = X_ED(iClass{k},:);
                myNS{iq,it}.mesh{k}.Y_ED = Y_ED(iClass{k},iq,it);
                myNS{iq,it}.mesh{k}.eLOO = inf;
                myNS{iq,it}.mesh{k}.class = k;
                myNS{iq,it}.mesh{k}.prob = probClass{k};
            end
        end

    % Otherwise, initialize single element mesh
    else        
        for iq = 1:nq
            myNS{iq,it}.mesh{1}.a = P_lb;
            myNS{iq,it}.mesh{1}.b = P_ub;
            myNS{iq,it}.mesh{1}.N_ED = size(X_ED,1);
            myNS{iq,it}.mesh{1}.X_ED = X_ED;
            myNS{iq,it}.mesh{1}.Y_ED = Y_ED(:,iq,it);
            myNS{iq,it}.mesh{1}.eLOO = inf;
            myNS{iq,it}.mesh{1}.class = [];
            myNS{iq,it}.mesh{1}.prob = 1;
        end
        fprintf('\nUsing single element.\n')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%% Adaptive pce %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Loop over quantities of interest
    for iq = 1:nq
        % Loop over elements
        for k = 1:myNS{iq,it}.Ne
            % Record time to fit pce model in element k
            tfit_k_start = tic;

            % Create local input object in uqlab
            for i = 1:nparam
                Input.Marginals(i).Type = param_distribution;
                Input.Marginals(i).Parameters = [myNS{iq,it}.mesh{k}.a(i), myNS{iq,it}.mesh{k}.b(i)];
            end
            myInput = uq_createInput(Input);
            myNS{iq,it}.mesh{k}.myInput = myInput;
            metaopts.Input = myInput;
            
            % Run until error tolerance is met
            while myNS{iq,it}.mesh{k}.eLOO > eLOO_tol
                % Add samples to ED after first iteration or ED too small
                if myNS{iq,it}.mesh{k}.eLOO ~= inf || myNS{iq,it}.mesh{k}.N_ED < N_add
                    % Sample full space for global
                    if isempty(myNS{iq,it}.mesh{k}.class)
                        X_ED_new = uq_getSample(myNS{iq,it}.mesh{k}.myInput, N_add, sample_type); 
                    % Sample within class for ns
                    else
                        X_ED_new = getSamplesClass(N_add, myNS{iq,it}.td, tlist(it), myNS{iq,it}.mesh{k}.class, sample_type);
                    end
                    % Integrate the full dfba model at new samples
                    [S_ED_new,T_ED_new] = evalForwardModel(X_ED_new);
                    % Evaluate the qoi for all time points and all new samples
                    Y_ED_new = zeros(N_add,nq,ntime);
                    for it2 = 1:ntime
                        Y_ED_new(:,:,it2) = qoi(S_ED_new,T_ED_new,tlist(it2));
                    end
                    % Add samples to global ED
                    N_ED = N_ED + N_add;
                    X_ED = [X_ED ; X_ED_new];
                    Y_ED = [Y_ED ; Y_ED_new];
                    % Add samples to local ED
                    if isempty(myNS{iq,it}.mesh{k}.class)
                        for iq2 = 1:nq
                            myNS{iq2,it}.mesh{k}.N_ED = myNS{iq2,it}.mesh{k}.N_ED + N_add;
                            myNS{iq2,it}.mesh{k}.X_ED = [myNS{iq2,it}.mesh{k}.X_ED ; X_ED_new];
                            myNS{iq2,it}.mesh{k}.Y_ED = [myNS{iq2,it}.mesh{k}.Y_ED ; Y_ED_new(:,iq2,it)];
                        end
                    else
                        % Evaluate the discontinuity time for all new samples and store them
                        Td_ED_new = tdisc(S_ED_new,T_ED_new);
                        Td_ED = [Td_ED ; Td_ED_new];
                        % Classify the new ED points
                        iClassNew{1} = find(Td_ED_new(:,myNS{iq,it}.idisc) > tlist(it));
                        iClassNew{2} = find(Td_ED_new(:,myNS{iq,it}.idisc) <= tlist(it));
                        % Store in local ED
                        for k2 = 1:myNS{iq,it}.Ne
                            for iq2 = 1:nq
                                myNS{iq2,it}.mesh{k2}.N_ED = myNS{iq2,it}.mesh{k2}.N_ED + length(iClassNew{k2});
                                myNS{iq2,it}.mesh{k2}.X_ED = [myNS{iq2,it}.mesh{k2}.X_ED ; X_ED_new(iClassNew{k2},:)];
                                myNS{iq2,it}.mesh{k2}.Y_ED = [myNS{iq2,it}.mesh{k2}.Y_ED ; Y_ED_new(iClassNew{k2},iq2,it)];
                            end
                        end
                    end
                end
                
                % If maximum deviation is less than tolerance, use constant pce
                if sqrt(var(myNS{iq,it}.mesh{k}.Y_ED)) <= std_factor*eLOO_tol
                    % Fit with a constant PCE
                    metaopts_cons = metaopts;
                    metaopts_cons.Degree = 1;
                    metaopts_cons.TruncOptions.qNorm = 1.0;
                    metaopts_cons.ExpDesign.X = myNS{iq,it}.mesh{k}.X_ED;
                    metaopts_cons.ExpDesign.Y = repmat(mean(myNS{iq,it}.mesh{k}.Y_ED),[myNS{iq,it}.mesh{k}.N_ED,1]);
                    myNS{iq,it}.mesh{k}.myPCE = uq_createModel(metaopts_cons);                    
                    % Set the LOO error to be zero
                    myNS{iq,it}.mesh{k}.myPCE.Error.ModifiedLOO = 0;                    
                    % Print statement
                    fprintf('\nConstant PCE used with LOO error estimate: 0\n')
                    
                % Otherwise, use the local ED
                else
                    % Create and add to uqlab the PCE model
                    metaopts.ExpDesign.X = myNS{iq,it}.mesh{k}.X_ED;
                    metaopts.ExpDesign.Y = myNS{iq,it}.mesh{k}.Y_ED;
                    myNS{iq,it}.mesh{k}.myPCE = uq_createModel(metaopts);
                end
                
                % Extract the modified LOO error
                myNS{iq,it}.mesh{k}.eLOO = myNS{iq,it}.mesh{k}.myPCE.Error.ModifiedLOO;
                
                % Break out of loop if max number of samples exceeded
                if myNS{iq,it}.mesh{k}.N_ED >= floor(N_max*myNS{iq,it}.mesh{k}.prob)
                    break
                end               
            end
            
            % Record if tolerance is met within element k
            if myNS{iq,it}.mesh{k}.eLOO <= eLOO_tol
                myNS{iq,it}.mesh{k}.meetTol = 1;
            else
                myNS{iq,it}.mesh{k}.meetTol = 0;
            end
            
            % Record time to fit pce in element k
            myNS{iq,it}.mesh{k}.tfit = toc(tfit_k_start);
            
            % Print status
            fprintf('\n')
            fprintf('finished element %g / %g for qoi %g / %g and time %g / %g\n',k,myNS{iq,it}.Ne,iq,nq,it,ntime)
            fprintf('total ED so far = %g\n',N_ED)
            fprintf('\n')
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



%% Post-processing of model

% Update Td models with all available data
if strcmp(fit_type, 'ns')
    for i = 1:ndisc
        metaopts.Input = tdisc_pce{i}.myInput;
        metaopts.ExpDesign.X = X_ED;
        if size(X_ED,1) == size(Td_ED(:,i),1)
            metaopts.ExpDesign.Y = Td_ED(:,i);
            tdisc_pce{i}.myPCE = uq_createModel(metaopts);
            tdisc_pce{i}.eLOO = tdisc_pce{i}.myPCE.Error.ModifiedLOO;
        end
    end
    for it = 1:ntime
        for iq = 1:nq
            if strcmp(myNS{iq,it}.fit_type, 'ns')
                myNS{iq,it}.td = tdisc_pce{myNS{iq,it}.idisc};
            end
        end
    end
end

% Compare surrogate to true model at validation points
if ~isempty(validation_mat_file)
    % Load the validation data
    load(validation_mat_file, 'X_val', 'Y_val')
    
    % Initialize variables
    N = size(X_val,1);
    Y_pce = zeros(N,nq,ntime);
    RMSE_val = zeros(nq,ntime);
    
    % Loop over time points
    for it = 1:ntime
        % Loop over qoi's
        for iq = 1:nq
            % Classify validation points in regions
            if strcmp(myNS{iq,it}.fit_type, 'ns')
                Td_val = uq_evalModel(myNS{iq,it}.td.myPCE, X_val);
                indexRegion{1} = find(Td_val > tlist(it));
                indexRegion{2} = find(Td_val <= tlist(it));                
            else
                indexRegion{1} = 1:N;
            end
            % Evaluate the surrogate
            for k = 1:myNS{iq,it}.Ne
                Y_pce(indexRegion{k},iq,it) = uq_evalModel(myNS{iq,it}.mesh{k}.myPCE, X_val(indexRegion{k},:));
            end
            % Enforce positivity
            index_neg = find(Y_pce(:,iq,it) < 0);
            Y_pce(index_neg,iq,it) = zeros(length(index_neg),1);
            % Calculate relative mean square error (RMSE)
            RMSE_val(iq,it) = mean((Y_val(:,iq,it)-Y_pce(:,iq,it)).^2) / var(Y_val(:,iq,it),1);
        end
    end
    
    % Plot parity
    for iq = 1:nq
        figure; hold on;
        for it = 1:ntime
            plot(Y_val(:,iq,it), Y_pce(:,iq,it), 'b.', 'MarkerSize', 15)
            plot([min(Y_val(:,iq,it)), max(Y_val(:,iq,it))],[min(Y_val(:,iq,it)), max(Y_val(:,iq,it))],'-k')
            set(gcf,'color','w')
            set(gca,'FontSize',20)
            xlabel('exact model')
            ylabel('surrogate model')
        end
        title(['qoi ' num2str(iq)])
    end
end

% Save model
if ~isempty(save_mat_file)
    save(save_mat_file)
end

% Remove added folders from path 
rmpath('./dfba_model')
