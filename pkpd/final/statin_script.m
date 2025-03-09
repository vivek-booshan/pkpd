% Parameter values for ihatemylifeODE
p = struct();

% Volume terms (assumed in liters)
p.V = struct('GI', 1, 'peripheral', 10, 'liver', 1);

% GI parameters (absorption rates, assumed in 1/hr)
p.GI = struct('chylomicron', struct('chol', 0.1, 'TG', 0.1));

% Peripheral tissue parameters (transport and clearance rates, assumed in 1/hr)
p.peripheral = struct( ...
    'chylomicron', struct('chol', 0.05, 'TG', 0.05), ...
    'HDL', struct('chol', 0.03, 'TG', 0.02), ...
    'clearance', struct('chol', 0.01, 'TG', 0.01), ...
    'LDL', struct('chol', 0.07) ... % Peripheral LDL cholesterol transport rate (1/hr)
);

% Liver parameters (transport and clearance rates, assumed in 1/hr)
p.liver = struct( ...
    'LDL', struct('chol', 0.07, 'TG', 0.05), ...
    'HDL', struct('chol', 0.04, 'TG', 0.03), ...
    'clearance', struct('chol', 0.02, 'TG', 0.01, 'statin', 0.05, 'ezetimibe', 0.05), ...
    'oxLDL', struct('chol', 0.005) ...
);

% ROS dynamics parameters
p.basalROS = 0.02;                % Baseline ROS production (arbitrary units/hr)
p.foodProductionMultiplier = 1.2; % ROS increase due to food intake (unitless multiplier)
p.antioxidant = 0.01;             % Antioxidant clearance of ROS (1/hr)
p.foodAntiOxMultiplier = 1.0;     % Antioxidant effect (unitless multiplier)

% Statin parameters
p.F = 0.15;       % Bioavailability of statin (arbitrary assumption)
p.ka = 1.0;       % Absorption rate constant for statin (1/hr)
p.liver_statin = 0.1;   % Liver statin concentration (arbitrary units)
p.liver_ezetimibe = 0.1;

% Initial conditions (assumed in arbitrary units or mg/dL)
y0 = [
    100;    % GI_chol
    50;     % Peripheral_chol
    20;     % Liver_chol
    0;      % Clearance_chol
    200;    % GI_TG
    100;    % Peripheral_TG
    50;     % Liver_TG
    0;      % Clearance_TG
    1;      % ROS (arbitrary baseline)
    0;      % oxLDL (starts at 0)
    10;     % Oral statin (arbitrary starting amount)
    0;      % Statin (arbitrary starting amount)
    0;    % Liver statin (arbitrary starting amount)
    0;      % Clearance statin (starts at 0)
]';
%%

% Run simulations
tspan = 1:1/60:24; % Simulate for 24 hours with 1-minute resolution
[~, y_nodrug] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0(1:10)); % No drug model
[~, y_statin] = ode45(@(t, y) statinODE(t, y, p), tspan, y0);       % Statin model
[t, y_ezetimibe] = ode45(@(t, y) ezetimibeODE(t, y, p), tspan, y0); % Ezetimibe model

% Extract relevant indices for cholesterol and oxLDL
chol_indices = 1:4; % GI_chol, Peripheral_chol, Liver_chol, Clearance_chol
TG_indices = 5:8; % GI_TG, Peripheral_TG, Liver_TG, Clearance_TG
oxLDL_index = 10;            % oxLDL

% Plotting setup
figure;
tiledlayout(2, 1); % Create two rows of plots: one for cholesterol, one for oxLDL

% Cholesterol dynamics
nexttile;
hold on;
plot(t, y_nodrug(:, chol_indices), 'LineWidth', 1.5); % No drug model
plot(t, y_statin(:, chol_indices), '--', 'LineWidth', 1.5); % Statin model
plot(t, y_ezetimibe(:, chol_indices), ':', 'LineWidth', 1.5); % Ezetimibe model
xlabel('Time (hours)');
ylabel('Cholesterol (mg/dL)');
title('Cholesterol Dynamics');
legend({'GI (No Drug)', 'Peripheral (No Drug)', 'Liver (No Drug)', 'Clearance (No Drug)', ...
        'GI (Statin)', 'Peripheral (Statin)', 'Liver (Statin)', 'Clearance (Statin)', ...
        'GI (Ezetimibe)', 'Peripheral (Ezetimibe)', 'Liver (Ezetimibe)', 'Clearance (Ezetimibe)'}, ...
        'Location', 'northeastoutside');
grid on;

% oxLDL dynamics
nexttile;
hold on;
plot(t, y_nodrug(:, oxLDL_index), 'LineWidth', 1.5); % No drug model
plot(t, y_statin(:, oxLDL_index), '--', 'LineWidth', 1.5); % Statin model
plot(t, y_ezetimibe(:, oxLDL_index), ':', 'LineWidth', 1.5); % Ezetimibe model
xlabel('Time (hours)');
ylabel('oxLDL (arbitrary units)');
title('Oxidized LDL Dynamics');
legend({'No Drug', 'Statin', 'Ezetimibe'}, 'Location', 'best');
grid on;

% Enhance overall figure aesthetics
sgtitle('Comparison of Cholesterol and oxLDL Dynamics');

%%
% Define color vector for cholesterol compartments
chol_colors = lines(4); % Single color from MATLAB's "lines" colormap

% Define line styles for each treatment model
line_styles = {'-', '--', ':'}; % No drug, statin, ezetimibe

% Plotting setup
figure;
tiledlayout(2, 1); % Create two rows of plots: one for cholesterol, one for oxLDL

% Cholesterol dynamics
nexttile;
hold on;

% Plot for each treatment model
models = {y_nodrug, y_statin, y_ezetimibe};
for i = 1:3
    for j = 1:length(TG_indices) % Iterate over compartments (GI, peripheral, liver, clearance)
        plot(t, models{i}(:, TG_indices(j)), ...
            'LineStyle', line_styles{i}, ...
            'Color', chol_colors(j, :), ...
            'LineWidth', 1.5);
    end
end

xlabel('Time (hours)');
ylabel('TG (mg/dL)');
title('Triglycerides Dynamics');
legend({'GI (No Drug)', 'Peripheral (No Drug)', 'Liver (No Drug)', 'Clearance (No Drug)', ...
        'GI (Statin)', 'Peripheral (Statin)', 'Liver (Statin)', 'Clearance (Statin)', ...
        'GI (Ezetimibe)', 'Peripheral (Ezetimibe)', 'Liver (Ezetimibe)', 'Clearance (Ezetimibe)'}, ...
        'Location', 'eastoutside');
grid on; 

% oxLDL dynamics
nexttile;
hold on;

% Plot oxLDL dynamics for each model
for i = 2:3
    plot(t, models{i}(:, oxLDL_index), ...
         'LineStyle', line_styles{i}, ...
        'Color', [0, 0.5, 0], ... % Example: green color for oxLDL
        'LineWidth', 1.5);
end

xlabel('Time (hours)');
ylabel('oxLDL (arbitrary units)');
title('Oxidized LDL Dynamics');
legend({'Statin', 'Ezetimibe'}, 'Location', 'best');
grid on; 

% Enhance overall figure aesthetics
sgtitle('Comparison of Triglyceride and oxLDL Dynamics');

%%
% Define indices for the statin variables
statin_indices = [11, 12, 13, 14]; % Oral statin, systemic statin, liver statin, clearance statin

% Define colors for the compartments
drug_colors = lines(length(statin_indices)); % Generate 4 distinct colors

% Create figure
figure;
hold on;

% Plot statin concentrations
for i = 1:length(statin_indices)
    plot(t, y_statin(:, statin_indices(i)), ...
        'LineWidth', 1.5, ...        % Line thickness
        'Color', drug_colors(i, :)); % Assign unique color to each compartment
end

% Add labels and legend
xlabel('Time (hours)');
ylabel('Drug Concentration (mg/L)');
title('Drug Concentrations');

% Generate legend dynamically
legend_entries = ["Oral", "Systemic", "Liver", "Clearance"];
legend(legend_entries, 'Location', 'northeast');

grid on;
hold off;
