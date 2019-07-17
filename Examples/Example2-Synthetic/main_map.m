
%% Initialize

% Start uqlab -- NOTE: MUST DOWNLOAD UQLAB (free for academic users https://www.uqlab.com/)
uqlab

% Clear the workspace
clear

% Fix the random seed to get repeatable results
rng(120)



%% User inputs

% Load the model
load('ns_pce_tol=1e-3', 'myNS', 'tlist', 'P_lb', 'P_ub', 'param_distribution', 'qoi', 'nq')

% Add dfba integrator to path -- NOTE: MUST BE MODIFIED BY USER (see evalForwardModel.m for more information)
addpath('./dfba_model')

% Add solvopt to path -- NOTE: MUST DOWNLOAD SOLVOPT (https://imsc.uni-graz.at/kuntsevich/solvopt/download.html#Matlab)
addpath('/Users/Joel/Documents/Matlab_toolboxes/so11matl')

% Specify parameter prior properties
x0 = (P_lb+P_ub)/2;
x_factor = 0.1;
Cx = diag((x0*x_factor).^2);

% Specify the noise properties
Cy = diag(([0.02, 0.4, 0.01, 0.05, 0.01, 0.1, 0.1]/2).^2);

% Specify if full model (1) or surrogate model (0)
full_model = 0;
if full_model == 1
    model = @(x)evalForwardModelQOI(x, tlist);
else
    model = @(x)evalSurrogateModelFull(x, myNS, tlist);
end

% Set solver options
solvopt_options = SOPTIONS;
solvopt_options(3) = 1e-4;
solvopt_options(5) = 1;

% Save the workspace? Empty vector = No
save_mat_file = [];



%% Get simulated data

% List of time points
T = tlist';
nT = length(T);

% Specify marginals
nparam = length(P_lb);
for i = 1:nparam
    Input.Marginals(i).Type = 'Gaussian';
    Input.Marginals(i).Parameters = [x0(i), x0(i)*x_factor];
end

% Create and add the resulting input object to uqlab
myInput = uq_createInput(Input);

% Specify true parameter value
x_true = [1.55e0, 4.7e-2, 1.53e1, 2.65e-1, 4.64e-1, 2.15e0, 1.28e0, 1.9e-1, 0.94e-2, 1.61e1, 2.78e-1, 0.95e0, 3.7e0, 5.35e-1, 1.54e0, 1.05e0, 1.05e0, 2.15e0, 1.0e0, 2.0e0];

% Calculate outputs at true parameter value
[S_true,T_true] = evalForwardModel(x_true);
y_true = zeros(nT,nq);
for i = 1:nT
    y_true(i,:) = qoi(S_true,T_true,T(i));
end

% Get noise realization and add to true output to get noisy data
error = zeros(nT,nq);
for i = 1:nT
    for j = 1:nq
        error(i,j) = (Cy(j,j))*randn;
    end
end
y = y_true + error;

% Check if negative, and set them to zero
for i = 1:nT
    for j = 1:nq
        if y(i,j) < 0
            y(i,j) = 0.0;
        end
    end
end

% Change dimensions
error = error';
y = y';



%% Run optimization

% Get dimensions
nx = length(x_true);
ny = size(y,1);
N = size(y,2);

% Specify the map objective function
fun = @(x)(map_objective_func(x, y, inv(Cy), x0, inv(Cx), model, zeros(ny,N), zeros(1,nx)));

% Define the gradient function using uqlab if surrogate
if full_model == 0
    grad = @(x)(uq_gradient(x, fun, 'forward', 'relative', 1e-4, myInput.Marginals));
    hessian = @(x)(uq_gradient(x, grad, 'forward', 'relative', 1e-4, myInput.Marginals));
end

% % Constraints on bounds
% func = @(x)(calc_max_residual_bounds(x, P_lb, P_ub));

% Use solvopt (with finite differences) to solve MAP problem
tic
if full_model == 1
    [x,f,solvopt_options] = SOLVOPT(x0, fun, [], solvopt_options, [], []);
else
    [x,f,solvopt_options] = SOLVOPT(x0, fun, grad, solvopt_options, [], []);    
end
x = abs(x);
toc

% Get model gradient
if full_model == 1
    G = zeros(ny*N,nx);
    m = @(x)(reshape(model(x),[ny*N,size(x,1)]));
    epsVal = 1e-4;
    f = m(x);
    for i = 1:nx
        x_new = x;
        x_new(i) = x(i) + epsVal*x(i);
        f_new = m(x_new);
        G(:,i) = (f_new - f)/(epsVal*x(i));
    end
else
    m = @(x)(reshape(model(x),[size(myNS,1)*size(myNS,2),size(x,1)])');
    G = uq_gradient(x, m, 'forward', 'relative', 1e-4, myInput.Marginals);
    G = squeeze(G)';
end

% Calculate the covariance of the parameter estimates
p_map = x';
Cy_full = Cx;
G = [eye(nx) ; G];
Cy_full = blkdiag(Cy_full,Cy);
for i = 1:N-1
    Cy_full = blkdiag(Cy_full, Cy);
end
Cov = inv(G'*inv(Cy_full)*G);
cov_diag = diag(Cov);
p_std = sqrt(cov_diag);

% Simulate system under map
[S0,T0] = evalForwardModel(x0);
[S_map,T_map] = evalForwardModel(p_map');



%% Generate plots

% Plot the posterior
figure('units','normalized','outerposition',[0 0 1 1])
%set the scale factor around 1
scale = 0.3;
%choose colormap
CM = brighten(jet(nx),0);
CM(13,:) = [1, 0.9, 0];
for i = 1:nx
    % Select subplot
    subplot(nx, nx, i+nx*(i-1))
    hold on
    
    % Calculate posterior 
    mean_i = p_map(i)/x0(i);
    std_i = sqrt(Cov(i,i)/x0(i)^2);
    x_list = linspace(mean_i-4*std_i, mean_i+4*std_i, 100); 
    y_list = normpdf(x_list, mean_i, std_i);
    plot(x_list, y_list, '-', 'color', CM(i,:), 'LineWidth', 2)
        
    % Set the axis
    axis([1-scale, 1+scale, 0, 1.5*max(y_list)])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    % Loop over parameters up to current one
    for j = 1:i-1
        % Select subplot
        subplot(nx, nx, j+nx*(i-1))
        hold on
        
        % Calculate the mean and covariance between pi, pj
        mean_ji = [p_map(j)/x0(j) ; p_map(i)/x0(i)];
        Cov_ji = Cov([j;i],[j;i]);
        Cov_ji(1,1) = Cov_ji(1,1)/x0(j)^2;
        Cov_ji(2,2) = Cov_ji(2,2)/x0(i)^2;
        Cov_ji(1,2) = Cov_ji(1,2)/x0(j)/x0(i);
        Cov_ji(2,1) = Cov_ji(2,1)/x0(j)/x0(i);
        plotErrorEllipse(mean_ji, Cov_ji, 0.95, CM(i,:), 2)
        
        % Plot the true value
        plot(x_true(j)/x0(j), x_true(i)/x0(i), 'kx', 'Markersize', 5, 'Linewidth', 2)
        
        % Set the axis
        axis([1-scale 1+scale 1-scale 1+scale])
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
end
%define labels for all parameters
p_labels{1} = 'v_{max,C}';
p_labels{2} = 'K_{C}';
p_labels{3} = 'K_{i,E}';
p_labels{4} = 'v_{max,N}';
p_labels{5} = 'K_{N}';
p_labels{6} = 'v_{max,O}';
p_labels{7} = 'K_{O}';
p_labels{8} = 'v_{ATP,m}';
p_labels{9} = 'X_0';
p_labels{10} = 'C_0';
p_labels{11} = 'N_0';
p_labels{12} = 'O_0';
p_labels{13} = 'S_{C,vx}';
p_labels{14} = 'S_{N,vx}';
p_labels{15} = 'S_{ATP,vx}';
p_labels{16} = 'S_{ATP,vox}';
p_labels{17} = 'S_{ATP,vf}';
p_labels{18} = 'S_{ATP,vlip}';
p_labels{19} = 'S_{OX,vox}';
p_labels{20} = 'S_{OX,vf}';
%put labels along first column and bottom row
for i = 1:nx
    subplot(nx, nx, 1+nx*(i-1))
    ylabel(p_labels{i}, 'Color', CM(i,:), 'fontweight','bold')
    ylh = get(gca,'ylabel');
    ylp = get(ylh, 'Position');
    ext = get(ylh,'Extent');
    set(ylh, 'Rotation',0, 'Position',ylp-[ext(3) 0 0], 'VerticalAlignment','middle')
    subplot(nx, nx, i+nx*(nx-1))
    xlabel(p_labels{i}, 'fontweight','bold')
end
set(gcf,'color','w')

% Plot the simulated data (biomass/lipids) with predicted dfba state profiles
Y = y';
figure
hold on
plot(T_map{1},S_map{1}(:,1),'-b','LineWidth',2)
plot(T_map{1},S_map{1}(:,5),'-r','LineWidth',2)
plot(T0{1},S0{1}(:,1),'--b','LineWidth',2)
plot(T0{1},S0{1}(:,5),'--r','LineWidth',2)
plot(T,Y(:,1),'bs','MarkerSize',8,'MarkerFaceColor','b')
plot(T,Y(:,5),'rs','MarkerSize',8,'MarkerFaceColor','r')
set(gcf,'color','w')
set(gca,'FontSize',16)
xlabel('time (hr)')
ylabel('concentration (g/L)')
legend('biomass', 'lipids')
axis([0 40 0 0.7])

% Plot the simulated data (species) with predicted dfba state profiles
Y = y';
figure
hold on
plot(T_map{1},S_map{1}(:,2),'-b','LineWidth',2)
plot(T_map{1},S_map{1}(:,3),'-r','LineWidth',2)
plot(T_map{1},S_map{1}(:,4),'-g','LineWidth',2)
plot(T_map{1},S_map{1}(:,6),'-k','LineWidth',2)
plot(T_map{1},S_map{1}(:,7),'-m','LineWidth',2)
plot(T0{1},S0{1}(:,2),'--b','LineWidth',2)
plot(T0{1},S0{1}(:,3),'--r','LineWidth',2)
plot(T0{1},S0{1}(:,4),'--g','LineWidth',2)
plot(T0{1},S0{1}(:,6),'--k','LineWidth',2)
plot(T0{1},S0{1}(:,7),'--m','LineWidth',2)
plot(T,Y(:,2),'bs','MarkerSize',8,'MarkerFaceColor','b')
plot(T,Y(:,3),'rs','MarkerSize',8,'MarkerFaceColor','r')
plot(T,Y(:,4),'gs','MarkerSize',8,'MarkerFaceColor','g')
plot(T,Y(:,6),'ks','MarkerSize',8,'MarkerFaceColor','k')
plot(T,Y(:,7),'ms','MarkerSize',8,'MarkerFaceColor','m')
set(gcf,'color','w')
set(gca,'FontSize',16)
xlabel('time (hr)')
ylabel('concentration (g/L)')
legend('carbon', 'nitrogen', 'oxygen', 'ethanol', 'oxidation')
axis([-1 40 0 16.5])



%% Save the model

% Save model
if ~isempty(save_mat_file)
    save(save_mat_file)
end

% Remove added folders from path 
rmpath('./dfba_model')
rmpath('/Users/Joel/Documents/Matlab_toolboxes/so11matl')
