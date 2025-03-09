%% Cholesterol and oxLDL Dynamics with Niacin and PCSK9 Inhibitor

% Define color vector for cholesterol compartments
chol_colors = lines(4); % 4 compartments × 3 models

% Define line styles for each treatment model
line_styles = {'-', '--', ':'}; % No drug, Niacin, PCSK9 Inhibitor

% Initialize Parameters and Initial Conditions
tspan = 1:1/60:24; % Simulate for 24 hours with 1-minute resolution
y0 = [100, 80, 50, 10, 200, 150, 100, 50, 0, 1000]; % Example initial values
% Parameter Initialization
p = struct();

% Define parameter structure 'p'
p.V.GI = 1.0;                   % Volume of GI compartment
p.V.peripheral = 1.0;           % Volume of peripheral compartment
p.V.liver = 1.0;                % Volume of liver compartment
p.V.clearance = 1.0;            % Volume of clearance compartment

% GI parameters
p.GI.chylomicron.chol = 0.1;    % Rate constant for cholesterol
p.GI.chylomicron.TG = 0.2;      % Rate constant for triglycerides

% Liver parameters
p.liver.LDL.chol = 0.3;         % LDL cholesterol transport rate
p.liver.HDL.chol = 0.4;         % HDL cholesterol transport rate
p.liver.clearance.chol = 0.02;   % Liver cholesterol clearance rate
p.liver.oxLDL.chol = 0.1;       % oxLDL generation rate from liver cholesterol
p.liver.LDL.TG = 0.2;           % LDL triglyceride transport rate
p.liver.HDL.TG = 0.3;           % HDL triglyceride transport rate
p.liver.clearance.TG = 0.01;     % Liver triglyceride clearance rate

% Peripheral parameters
p.peripheral.LDL.chol = 0.3;
p.peripheral.HDL.chol = 0.4;
p.peripheral.clearance.chol = 0.01;
p.peripheral.chylomicron.chol = 0.2;
p.peripheral.chylomicron.TG = 0.2;
p.peripheral.HDL.TG = 0.1;
p.peripheral.clearance.TG = 0.01;

% ROS dynamics parameters
p.basalROS = 0.01;              % Basal ROS production rate
p.foodProductionMultiplier = 1; % Multiplier for food-induced ROS
p.antioxidant = 0.005;          % Antioxidant effect rate
p.foodAntiOxMultiplier = 1.5;   % Multiplier for food-induced antioxidant effect


% Solve ODEs for each treatment model
[t_nodrug, y_nodrug] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0); % No drug model
[t_niacin, y_niacin] = ode45(@(t, y) niacinODE(t, y, p), tspan, y0);
[t_pcsk9, y_pcsk9] = ode45(@(t, y) pcsk9ODE(t, y, p), tspan, y0);

% Combine models for plotting
models = {y_nodrug, y_niacin, y_pcsk9};
t = t_nodrug; % Assuming all models use the same time span

%%
% Extract relevant indices
chol_indices = [1:4];  % Cholesterol compartments: GI, Peripheral, Liver, Clearance
TG_indices = 5:8;    % TG compartments
oxLDL_index = 10;    % oxLDL
chol_colors = lines(5);

% Plotting setup
figure(1); clf
% tiledlayout(2, 1); % Create two rows of plots: one for cholesterol, one for oxLDL

% nexttile;
hold on;

% Plot for each treatment model
for i = 1:3
    for j = 1:length(chol_indices) % Iterate over compartments (GI, peripheral, liver, clearance)
        plot(t, models{i}(:, chol_indices(j)), ...
            'LineStyle', line_styles{i}, ...
            'Color', chol_colors(j, :), ...
            'LineWidth', 1.5);
    end
end
xlabel('Time (hours)');
ylabel('Triglycerides (mg/dL)');
title('Triglyceride Dynamics');
% legendEntries = {};
% treatments = {'No Drug', 'Niacin', 'PCSK9'};
% compartments = {'GI', 'Peripheral', 'Liver', 'Clearance'};
legend({'GI (No Drug)', 'Peripheral (No Drug)', 'Liver (No Drug)', 'Clearance (No Drug)', ...
        'GI (Niacin)', 'Peripheral (Niacin)', 'Liver (Niacin)', 'Clearance (Niacin)', ...
        'GI (PCSK9)', 'Peripheral (PCSK9)', 'Liver (PCSK9)', 'Clearance (PCSK9)'}, ...
        'Location', 'eastoutside');
grid on;

% nexttile;
% hold on;
% 
% % Plot oxLDL dynamics for Niacin and PCSK9 models
% for i = 2:3 % Niacin and PCSK9 Inhibitor
%     plot(t, models{i}(:, oxLDL_index), ...
%          'LineStyle', line_styles{i}, ...
%          'Color', [0, 0.5, 0], ... % Example: green color for oxLDL
%          'LineWidth', 1.5);
% end

% xlabel('Time (hours)');
% ylabel('oxLDL (arbitrary units)');
% title('Oxidized LDL Dynamics');
% legend({'Niacin', 'PCSK9 Inhibitor'}, 'Location', 'best');
% grid on;
%%
% Time span for simulation
% tspan = [0 100]; % Adjust time span as needed

% Initial conditions
% y0 = [1; 1; 1; 0; 1; 1; 1; 0; 1; 1]; % Example initial values

% Parameters (example values, adjust as needed)
p.V.GI = 1;
p.GI.chylomicron.chol = 0.1;
p.GI.chylomicron.TG = 0.1;
p.V.peripheral = 1;
p.peripheral.LDL.chol = 0.1;
p.peripheral.chylomicron.chol = 0.1;
p.peripheral.HDL.chol = 0.1;
p.peripheral.clearance.chol = 0.1;
p.peripheral.clearance.TG = 0.1;
p.V.liver = 1;
p.liver.LDL.chol = 0.1;
p.liver.HDL.chol = 0.1;
p.liver.clearance.chol = 0.1;
p.liver.clearance.TG = 0.1;
p.liver.oxLDL.chol = 0.1;
p.V.clearance = 1;
p.basalROS = 0.1;
p.foodProductionMultiplier = 1;
p.antioxidant = 0.1;
p.foodAntiOxMultiplier = 1;

% Solve ODE for each condition
[t1, y1] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0);
[t2, y2] = ode45(@(t, y) pcsk9ODE(t, y, p), tspan, y0);
[t3, y3] = ode45(@(t, y) niacinODE(t, y, p), tspan, y0);

% Extract oxLDL dynamics
oxLDL_nodrug = y1(:, 10);
oxLDL_pcsk9 = y2(:, 10);
oxLDL_niacin = y3(:, 10);

% Plot results
figure;
plot(t1, oxLDL_nodrug, 'r-', 'LineWidth', 2); hold on;
plot(t2, oxLDL_pcsk9, 'b--', 'LineWidth', 2);
plot(t3, oxLDL_niacin, 'g-.', 'LineWidth', 2);
xlabel('Time');
ylabel('oxLDL Concentration');
title('oxLDL Dynamics under Different Conditions');
legend({'No Drugs', 'PCSK9 Inhibitor', 'Niacin'}, 'Location', 'Best');
grid on;
% Time span and initial conditions
tspan = [0 200]; % Adjust time span as needed
y0 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]; % Initial values for state variables

% Parameter structure (example parameters; adjust as needed)
p.V.GI = 1;
p.GI.chylomicron.chol = 0.1;
p.V.liver = 1;
p.liver.LDL.chol = 0.2;
p.liver.HDL.chol = 0.1;
p.V.peripheral = 1;
p.peripheral.LDL.chol = 0.1;
p.peripheral.clearance.chol = 0.05;
p.basalROS = 0.02;
p.foodProductionMultiplier = 1.0;
p.foodAntiOxMultiplier = 0.8;
p.antioxidant = 0.05;

% Solve ODEs for each condition
[t_noDrug, y_noDrug] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0);
[t_niacin, y_niacin] = ode45(@(t, y) niacinODE(t, y, p), tspan, y0);
[t_pcsk9, y_pcsk9] = ode45(@(t, y) pcsk9ODE(t, y, p), tspan, y0);

% Extract oxLDL dynamics
oxLDL_noDrug = y_noDrug(:, 10);
oxLDL_niacin = y_niacin(:, 10);
oxLDL_pcsk9 = y_pcsk9(:, 10);

% Plot
figure;
plot(t_noDrug, oxLDL_noDrug, 'k-', 'LineWidth', 1.5); hold on;
plot(t_niacin, oxLDL_niacin, 'b--', 'LineWidth', 1.5);
plot(t_pcsk9, oxLDL_pcsk9, 'r-.', 'LineWidth', 1.5);
xlabel('Time');
ylabel('oxLDL');
legend({'No Drug', 'Niacin', 'PCSK9 Inhibitor'}, 'Location', 'Best');
title('oxLDL Dynamics Under Different Treatments');
grid on;
