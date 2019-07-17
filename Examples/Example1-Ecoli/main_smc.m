
%% Initialize

% Start uqlab -- NOTE: MUST DOWNLOAD UQLAB (free for academic users https://www.uqlab.com/)
uqlab

% Clear the workspace
clear

% Fix the random seed to get repeatable results
rng(120)



%% User inputs

% Load the model
load('nsPCE_tol=1e-3', 'myNS', 'tlist', 'P_lb', 'P_ub', 'param_distribution', 'qoi', 'nq')

% Add dfba integrator to path -- NOTE: MUST BE MODIFIED BY USER (see evalForwardModel.m for more information)
addpath('./dfba_model')

% Specify the noise properties -- normal(0,std), std = std_factor*model + std_floor
std_factor = [0.05, 0.05, 0.05];
std_floor = [0.01, 0.01, 0.01];

% Specify the likelihood function
l_func_k = @(y,x,k)(likelihood_func(y, x, k, std_factor, std_floor, myNS, tlist));

% Number of particles in smc
N = 1e5;

% Save the workspace? Empty vector = No
save_mat_file = [];



%% Get simulated data

% List of time points
T = tlist';
nT = length(T);

% Specify marginals
nparam = length(P_lb);
for i = 1:nparam
    Input.Marginals(i).Type = param_distribution;
    Input.Marginals(i).Parameters = [P_lb(i), P_ub(i)];
end

% Create and add the resulting input object to uqlab
myInput = uq_createInput(Input);

% Get true parameter value
X_true = [10.216284801235412, 0.002507499646600, 5.819954774961766, 0.015151837822482, 0.004829248945044, 13.167435058221075]; 
%X_true = uq_getSample(myInput, 1, 'MC');

% Calculate outputs at true parameter value
[S_true,T_true] = evalForwardModel(X_true);
Y_true = zeros(nT,nq);
for i = 1:nT
    Y_true(i,:) = qoi(S_true,T_true,T(i));
end

% Get noise realization and add to true output to get noisy data
E = zeros(nT,nq);
for i = 1:nT
    for j = 1:nq
        E(i,j) = (std_factor(j)*Y_true(i,j)+std_floor(j))*randn;
    end
end
Y = Y_true + E;



%% Run the smc algorithm

% Get samples from prior
X = uq_getSample(myInput, N, 'MC');

% Define uniform weights
W = ones(N,1)/N;

% Initialize data storage
X_data = zeros(N,nparam,nT+1);
W_data = zeros(N,1,nT+1);
X_data(:,:,1) = X;
W_data(:,:,1) = W;

% Print statement
fprintf('%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('STARTING SMC ALGORITHM\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%\n\n')

% Loop over time
for k = 1:nT
    % Print starting statement
    fprintf('starting sample %g / %g...',k,nT)
    tic
    
    % Calculate likelihood for samples
    L = l_func_k(Y(k,:), X, k);

    % Reweighting -- recursively update weights and normalize
    W = W.*L;
    W = W./sum(W);
    
    % Resample if effective number of particles is less than 85% of original
    Neff = 1/sum(W.^2);
    if Neff < 0.85*N
        fprintf('need to resample...')
        [X,W,~] = systematicResample(X,W);
    end
    
    % Store data
    X_data(:,:,k+1) = X;
    W_data(:,:,k+1) = W;
    
    % Print time
    fprintf('took %g seconds\n',toc)
end
fprintf('\n')

% Get smc estimates as mean of posterior 
X_opt = mean(X,1);

% Simulate true system at smc estimates
[S_opt,T_opt] = evalForwardModel(X_opt);
Y_opt = zeros(nq,nT);
for i = 1:nT
    Y_opt(:,i) = qoi(S_opt,T_opt,T(i));
end



%% Calculate predictive distributions

% Print statement
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('GETTING PRIOR PREDICTIVE\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%\n\n')

% Prior distribution
X_prior = X_data(:,:,1);
Y_prior = zeros(N,nq,nT);
for i = 1:nT
    tic
    fprintf('starting time %g / %g...',i,nT)
    Y_prior(:,:,i) = evalSurrogateModel(X_prior, myNS, i, T);
    index = find(Y_prior(:,2,i) > 8);
    Y_prior(index,2,i) = 8;
    fprintf('took %g seconds\n',toc)
end
fprintf('\n')

% Print statement
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
fprintf('GETTING POSTERIOR PREDICTIVE\n')
fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n')

% Posterior distribution
X_post = X_data(:,:,end);
Y_post = zeros(N,nq,nT);
for i = 1:nT
    tic
    fprintf('starting time %g / %g...',i,nT)
    Y_post(:,:,i) = evalSurrogateModel(X_post, myNS, i, T);
    index = find(Y_post(:,2,i) > 8);
    Y_post(index,2,i) = 8;
    fprintf('took %g seconds\n',toc)
end
fprintf('\n')



%% Generate plots

% Plot histograms of parameter posterior on top of prior
X_lb = P_lb;
X_ub = P_ub;
figure
for i = 1:nparam
    subplot(nparam, nparam, i+nparam*(i-1))
    histogram(X_prior(:,i), 'BinLimits', [X_lb(i), X_ub(i)], 'BinWidth', (X_ub(i)-X_lb(i))/30, 'FaceColor', 'g', 'EdgeColor', 'none', 'Normalization', 'probability')
    hold on;
    numbin = histogram(X_post(:,i), 'BinLimits', [X_lb(i), X_ub(i)], 'BinWidth', (X_ub(i)-X_lb(i))/30, 'FaceColor', 'b', 'EdgeColor', 'none', 'Normalization', 'probability');
    axis([0.95*X_lb(i), 1.05*X_ub(i), 0, 1.5*max(numbin.Values)])
    line([X_true(i), X_true(i)], [0,1], 'LineWidth', 2, 'Color', 'r');
    if i == 1
        set(gca,'xtick',[])    
        set(gca,'ytick',[]) 
    end
    if i > 1 && i < nparam
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    if i == nparam
        set(gca,'ytick',[])
    end
    for j = 1:i-1
        subplot(nparam, nparam, j+nparam*(i-1))
        hold on;
        scatter(X_prior(:,j),X_prior(:,i),'g.')
        scatter(X_post(:,j),X_post(:,i),'b.')
        axis([0.95*X_lb(j), 1.05*X_ub(j), 0.95*X_lb(i), 1.05*X_ub(i)])
        if j == 1 && i < nparam
            set(gca,'xtick',[])
        end
        if j > 1 && i < nparam
            set(gca,'xtick',[])
            set(gca,'ytick',[])
        end
        if j > 1 && i == nparam
            set(gca,'ytick',[])
        end
    end
end
set(gcf,'color','w')

% Plot the simulated data with predicted dfba state profiles
figure
hold on
plot(T_opt{1},S_opt{1}(:,3),'-b','LineWidth',2)
plot(T_opt{1},S_opt{1}(:,4),'-r','LineWidth',2)
plot(T_opt{1},S_opt{1}(:,2),'-g','LineWidth',2)
plot(T,Y(:,1),'bx','MarkerSize',10,'LineWidth',2)
plot(T,Y(:,2),'rx','MarkerSize',10,'LineWidth',2)
plot(T,Y(:,3),'gx','MarkerSize',10,'LineWidth',2)
set(gcf,'color','w')
set(gca,'FontSize',16)
xlabel('time (hr)')
ylabel('concentration (g/L)')
legend('glucose', 'xylose', 'biomass')
axis([0 9.0 -0.02 16])

% Plot prior predictive
qoi_names = {'glucose', 'xylose', 'biomass'};
colors = {'b', 'r', 'g'};
binlimits = [0 15 ; 0 8 ; 0 12];
binwidth = [0.3, 0.2, 0.2];
ymax = [0.15, 0.2, 0.25];
figure
for i = 1:nT
    for j = 1:nq
        subplot(nT,nq,(i-1)*nq+j)
        histogram(Y_prior(:,j,i), 'BinLimits', binlimits(j,:), 'BinWidth', binwidth(j), 'FaceColor', colors{j}, 'EdgeColor', 'none', 'Normalization', 'probability')
        set(gcf,'color','w')
        set(gca,'FontSize',16)
        xlim(binlimits(j,:))
        ylim([0, ymax(j)])
        if i < nT
            set(gca,'xtick',[])
        else
            xlabel([qoi_names{j} ' (g/L)'])
        end
        set(gca,'ytick',[])
    end
end

% Save model
if ~isempty(save_mat_file)
    save(save_mat_file)
end

% Remove added folders from path 
rmpath('./dfba_model')
