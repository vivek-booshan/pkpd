% Adjusted parameters for high saturated fat intake
p = struct();

% Volume terms (assumed in liters)
p.V = struct('GI', 1, 'peripheral', 10, 'liver', 1);

% GI parameters (increased absorption rates, assumed in 1/hr)
p.GI = struct('chylomicron', struct('chol', 0.5, 'TG', 0.5)); % High absorption

% Peripheral tissue parameters (slightly reduced clearance rates)
p.peripheral = struct( ...
    'chylomicron', struct('chol', 0.05, 'TG', 0.05), ...
    'LDL', struct('chol', 0.03), ...
    'HDL', struct('chol', 0.02, 'TG', 0.015), ...
    'clearance', struct('chol', 0.005, 'TG', 0.005) ... % Impaired clearance
);

% Liver parameters (reduced clearance rates, assumed in 1/hr)
p.liver = struct( ...
    'LDL', struct('chol', 0.1, 'TG', 0.08), ...
    'HDL', struct('chol', 0.03, 'TG', 0.02), ...
    'clearance', struct('chol', 0.01, 'TG', 0.005), ... % Reduced clearance
    'oxLDL', struct('chol', 0.01) ... % Increased oxLDL production
);

% ROS dynamics parameters (increased oxidative stress)
p.basalROS = 0.05;                % Increased baseline ROS production
p.foodProductionMultiplier = 2.0; % ROS increase due to high saturated fat intake
p.antioxidant = 0.005;            % Slightly impaired antioxidant activity
p.foodAntiOxMultiplier = 0.8;     % Reduced antioxidant effect

% Initial conditions (assumed in arbitrary units or mg/dL)
y0 = [
    200;    % GI_chol (higher due to saturated fat intake)
    30;    % Peripheral_chol
    30;     % Liver_chol
    0;      % Clearance_chol
    400;    % GI_TG (higher due to saturated fat intake)
    100;    % Peripheral_TG
    100;    % Liver_TG
    0;      % Clearance_TG
    1;      % ROS (arbitrary baseline)
    0       % oxLDL (starts at 0)
]';

% Parameters including glucose effects
p_fructose = p; 
p_fructose.fructose = struct( ...
    'liver_effect', 2, ...          % Glucose increases liver cholesterol synthesis
    'TG_effect', 1.2, ...
    'ROS_effect', 0.02, ...            % Glucose increases ROS production
    'clearance_effect', 0.02, ...
    'consumption_rate', 0.01 ...       % Glucose consumption rate
);

% Adjust initial conditions to include glucose
y0_fructose = [
    200;    % GI_chol (higher due to saturated fat intake)
    30;    % Peripheral_chol
    30;     % Liver_chol
    0;      % Clearance_chol
    400;    % GI_TG (higher due to saturated fat intake)
    100;    % Peripheral_TG
    100;    % Liver_TG
    0;      % Clearance_TG
    1;      % ROS (arbitrary baseline)
    0       % oxLDL (starts at 0)
    100     % fructose
]';
% Time span for simulation
tspan = 1:1/60:24; % 24 hours with 1-minute resolution
chol_indices = [1:4, 10];
chol_colors = lines(5);

[t_fructose, y_fructose] = ode45(@(t, y) nodrugfructoseODE(t, y, p_fructose), tspan, y0_fructose);

% Run simulation for low-fat intake
[t, y] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0);

% Plot results
figure(1); clf; hold on;
models = {y, y_fructose};
for i = 1:2
    for j = 1:length(chol_indices) % Iterate over compartments (GI, peripheral, liver, clearance)
        plot(t, models{i}(:, chol_indices(j)), ...
            'LineStyle', line_styles{i}, ...
            'Color', chol_colors(j, :), ...
            'LineWidth', 1.5);
    end
end
title('Comparison of Intake of Fructose vs No Fructose on Lipid Dynamics');
xlabel('Time (hours)');
ylabel('Concentration (mg/dL)');
legend(...
       'No Fru  - GI Chol', 'No Fru  - Peripheral Chol', 'No Fru  - Liver Chol', ...
       'No Fru  - Clearance Chol', 'No Fru Fat - oxLDL', ...
       'Fru  - GI Chol', 'Fru  - Peripheral Chol', 'Fru  - Liver Chol', ...
       'Fru  - Clearance Chol', 'Fru  - oxLDL' ...
       );
grid on;

