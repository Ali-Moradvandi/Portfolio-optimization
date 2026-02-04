%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) TU Delft All Rights Reserved
% Written by Ali Moradvandi
% For any correspondence: moradvandi@gmail.com

%% Introduction of code (purpose)
% This code is an example code for the portfolio optimization 
% for the transition of industrial chemical and petrochemical 
% clusters from fossil feedstocks for three plants namely Olefins,
% Methanol-topolefins, and CO2-toethylene.

% The details of the methodology can be found in XXX.

% This code is based on the data given for Scenario 1 (market-based economic
% paramaters) in paper XXX.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear

%% ECONOMIC PARAMETERS
% Number of processes
N = 3;           % N
% Specific CAPEX per plant 
TCC = [1031.93; 322.04; 646.92];
% Expected return on investment
ROI = [4.84; -2.13; -9.06];
% Std dev of ROI 
sigma = [1.77; 3.68; 1.45];
% Filling diagonal matrix of std devs
S = diag(sigma);
lower_vals = [ ...
   -0.28;     % row 2, col 1
    0.19;     % row 3, col 1
    0.11      % row 3, col 2
];
L = zeros(N);  
L(tril(true(N), -1)) = lower_vals;
Phi = L + L';       
Phi(1:N+1:end) = 1;     % correlation matrix (NxN) 

% Total investment budget
TI = 1.05*max(TCC);                        
% Lower bounds on capacity ratios
c_LL = [0; 0; 0];
% Upper bounds on capacity ratios
c_UL = [1; 1; 1]; 

Max_Cap  = [878.2; 314.4; 328.2];
req_production = 878.2; %642.6;

%% PARETO FRONT COMPUTATION
% Tradeoff parameter values
alpha_vals = linspace(0,1,20);        
% Optimal ROI and Risk over pareto
pareto_ROI = zeros(length(alpha_vals), 1); 
pareto_Risk = zeros(length(alpha_vals), 1);
% Optimal capacity ratios over pareto
pareto_c = zeros(length(alpha_vals), N);
% Optimal investment weights over pareto 
pareto_w = zeros(length(alpha_vals), N+1);  

% Optimization settings
options = optimoptions('fmincon', 'Display', 'off');

% Initial guess for capacities
c0 = ones(N,1) * TI / (sum(TCC));  % Evenly divide capital over plants
beta = 1;       % if not one means nonlinear scaling over capital expenditure

%% MAIN LOOP OVER ALPHA VALUES
for i = 1:length(alpha_vals)
    alpha = alpha_vals(i); 

    % Define objective function as a function of capacities
    obj = @(c) weighted_obj2(c, TCC, ROI, S, Phi, alpha, beta);

    % Inequality constraint: sum(c .* TCC) <= TI
    A = TCC';  % 1xN
    b = TI;

    % Equality constraint: sum(c .* Max_Cap) = req_production
    Aeq = [878.2; 314.4; 328.2]';
    beq = 878.2;

    % Solve optimization
    [c_opt, ~] = fmincon(obj, c0, A, b, Aeq, beq, c_LL, c_UL, [], options);

    % Compute derived weight vector w
    denom = sum(c_opt.^beta .* TCC);
    w_opt = (c_opt.^beta .* TCC) / denom;

    % Store results
    pareto_c(i,:) = c_opt';
    pareto_w(i,:) = [w_opt', denom];
    pareto_ROI(i) = ROI' * w_opt;
    pareto_Risk(i) = sqrt(w_opt' * S * Phi * S * w_opt);
end

%% PLOT PARETO FRONT
% assign number of pareto points to be plotted
q = 7;

f1 = figure(1);
hold on
plot(pareto_Risk, pareto_ROI, '-or', 'LineWidth', 2,'MarkerSize',5);
xlabel('Portfolio Risk','Interpreter','Latex','fontsize',9);
ylabel('Portfolio Return (%)','Interpreter','Latex','fontsize',9);
grid on;
f1.Position = [100 100 400 300];

f2 = figure(2);
b = bar(pareto_c(1:q,:), 'stacked');
threshold = 0.1; 
for i = 1:q
    cum_val = 0;
    for j = 1:N
        val = pareto_c(i, j);
        y = val / 2 + cum_val;
        if val >= threshold
            text(i, y, sprintf('$%.2f$', val), ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'latex', ...
                'FontSize', 7);
        end
        cum_val = cum_val + val;
    end
end
xticks([]);
xticklabels({});
light_green = [0.7 1 0.7]; 
light_green_2 = [0.4 1 0.4];
light_red = [1 0.5 0.5];
set(b(1), 'FaceColor', light_red);      % Red for Olefins
set(b(2), 'FaceColor', light_green);    % Light green for Methanol-to-olefins
set(b(3), 'FaceColor', light_green_2);  % Green for CO2-to-ethylene
ylabel('Plant Capacity ratio ($c_p$)','Interpreter','Latex','fontsize',9);
xlabel('Pareto Solution Index','Interpreter','Latex','fontsize',9);
legend ('Olefins', 'Methanol-to-olefins', 'CO2-to-ethylene', 'Interpreter','Latex','fontsize',8)
grid on;
f2.Position = [100 100 400 300];


f3 = figure(3);
b = bar(pareto_w(1:q,1:N), 'stacked');
threshold = 0.05; % 5% threshold
normalized_data = pareto_w(1:q, 1:N) ./ sum(pareto_w(1:q, 1:N), 2);
for i = 1:7
    cum_val = 0;
    for j = 1:N
        val = normalized_data(i,j);
        y = val / 2 + cum_val;
        if val >= threshold
            text(i, y, sprintf('$%.0f\\%%$', val*100), ...
                'HorizontalAlignment', 'center', ...
                'Interpreter', 'latex', ...
                'FontSize', 7);
        end
        cum_val = cum_val + val;
    end
end
xticks(1:q);
xticklabels(compose('%.1f', pareto_w(1:q,N+1)));
light_green = [0.7 1 0.7]; 
light_green_2 = [0.4 1 0.4];
light_red = [1 0.5 0.5];
light_red_2 = [1 0.3 0.3];
set(b(1), 'FaceColor', light_red); 
set(b(2), 'FaceColor', light_green); 
set(b(3), 'FaceColor', light_green_2); 
xticks(1:q);
xticklabels(compose('%.1f', pareto_w(1:7,N+1)));
ax = gca;
ax.XAxis.FontSize = 7;
ylabel('Investment Allocation ($w_p$)','Interpreter','Latex','fontsize',9);
xlabel('Total Investemet (MEUR) based on Pareto Solution Index','Interpreter','Latex','fontsize',9);
ylim([0 1.5])
legend ('Olefins', 'Methanol-to-olefins', 'CO2-to-ethylene', 'Interpreter','Latex','fontsize',8)
grid on;
f3.Position = [100 100 400 300];