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

% Parameters for low fatty acid intake
p_low = struct();

% Volume terms (assumed in liters, no change)
p_low.V = struct('GI', 1, 'peripheral', 10, 'liver', 1);

% GI parameters (reduced absorption rates, assumed in 1/hr)
p_low.GI = struct('chylomicron', struct('chol', 0.02, 'TG', 0.02)); % Lower absorption

% Peripheral tissue parameters (improved clearance rates, assumed in 1/hr)
p_low.peripheral = struct( ...
    'chylomicron', struct('chol', 0.05, 'TG', 0.05), ...
    'LDL', struct('chol', 0.01), ... % Slightly lower LDL transport
    'HDL', struct('chol', 0.04, 'TG', 0.03), ... % Higher HDL function
    'clearance', struct('chol', 0.02, 'TG', 0.02) ... % Improved clearance
);

% Liver parameters (improved clearance rates, assumed in 1/hr)
p_low.liver = struct( ...
    'LDL', struct('chol', 0.05, 'TG', 0.04), ...
    'HDL', struct('chol', 0.05, 'TG', 0.04), ... % Higher HDL efficiency
    'clearance', struct('chol', 0.01, 'TG', 0.02), ... % Improved clearance
    'oxLDL', struct('chol', 0.003) ... % Reduced oxLDL production
);

% ROS dynamics parameters (lower oxidative stress)
p_low.basalROS = 0.01;                % Lower baseline ROS production
p_low.foodProductionMultiplier = 1.0; % Minimal ROS increase due to food intake
p_low.antioxidant = 0.02;             % Enhanced antioxidant activity
p_low.foodAntiOxMultiplier = 1.2;     % Increased antioxidant effect

% Initial conditions (assumed in arbitrary units or mg/dL)
y0_low = [
    50;     % GI_chol (lower due to reduced fat intake)
    30;     % Peripheral_chol
    30;     % Liver_chol
    0;      % Clearance_chol
    100;    % GI_TG (lower due to reduced fat intake)
    100;     % Peripheral_TG
    25;     % Liver_TG
    0;      % Clearance_TG
    1;      % ROS (baseline similar to high-fat case)
    0       % oxLDL (starts at 0)
]';

%%
chol_indices = [1:4, 10];
chol_colors = lines(5); % Single color from MATLAB's "lines" colormap

% Define line styles for each treatment model
line_styles = {'-', '--', ':'}; % No drug, statin, ezetimibe


% Time span for simulation
tspan = 1:1/60:24; % 24 hours with 1-minute resolution

[t, y] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0);

% Run simulation for low-fat intake
[t_low, y_low] = ode45(@(t, y) nodrugODE(t, y, p_low), tspan, y0_low);

% Plot results
figure(1); clf; hold on;

models = {y, y_low};
for i = 1:2
    for j = 1:length(chol_indices) % Iterate over compartments (GI, peripheral, liver, clearance)
        plot(t, models{i}(:, chol_indices(j)), ...
            'LineStyle', line_styles{i}, ...
            'Color', chol_colors(j, :), ...
            'LineWidth', 1.5);
    end
end

% plot(t, y(:, [1:4, 10]), 'LineWidth', 1.5); hold on;
% plot(t_low, y_low(:, [1:4, 10]), '--', 'LineWidth', 1.5);
title('Comparison of High and Low Fat Intake on Lipid Dynamics');
xlabel('Time (hours)');
ylabel('Concentration (mg/dL)');
legend('High Fat - GI Chol', 'High Fat - Peripheral Chol', 'High Fat - Liver Chol', ...
       'High Fat - Clearance Chol', 'High Fat - oxLDL', ...
       'Low Fat - GI Chol', 'Low Fat - Peripheral Chol', 'Low Fat - Liver Chol', ...
       'Low Fat - Clearance Chol', 'Low Fat - oxLDL');
grid on;

