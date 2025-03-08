% Adjusted parameters for high saturated fat intake
p_PUFA = struct();

% Volume terms (assumed in liters)
p_PUFA.V = struct('GI', 1, 'peripheral', 10, 'liver', 1);

% GI parameters (increased absorption rates, assumed in 1/hr)
p_PUFA.GI = struct('chylomicron', struct('chol', 0.5, 'TG', 0.5)); % High absorption

% Peripheral tissue parameters (slightly reduced clearance rates)
p_PUFA.peripheral = struct( ...
    'chylomicron', struct('chol', 0.05, 'TG', 0.05), ...
    'LDL', struct('chol', 0.03), ...
    'HDL', struct('chol', 0.02, 'TG', 0.015), ...
    'clearance', struct('chol', 0.005, 'TG', 0.005) ... % Impaired clearance
);

% Liver parameters (reduced clearance rates, assumed in 1/hr)
p_PUFA.liver = struct( ...
    'LDL', struct('chol', 0.1, 'TG', 0.08), ...
    'HDL', struct('chol', 0.03, 'TG', 0.02), ...
    'clearance', struct('chol', 0.01, 'TG', 0.005), ... % Reduced clearance
    'oxLDL', struct('chol', 0.01) ... % Increased oxLDL production
);

p_PUFA.PUFA = struct(...
    'effectOnROS', -0.05, ...
    'effectOnHDL', 0.03, ...
    'effectOnLDL', -0.03, ...
    'metabolicClearance', 0.1, ...
    'giIntakeMultiplier', 1.5 ...
);

% ROS dynamics parameters (increased oxidative stress)
p_PUFA.basalROS = 0.05;                % Increased baseline ROS production
p_PUFA.foodProductionMultiplier = 2.0; % ROS increase due to high saturated fat intake
p_PUFA.antioxidant = 0.005;            % Slightly impaired antioxidant activity
p_PUFA.foodAntiOxMultiplier = 0.8;     % Reduced antioxidant effect

% Initial conditions (assumed in arbitrary units or mg/dL)
y0_PUFA = [
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

% Time span for simulation
tspan = 1:1/60:24; % 24 hours with 1-minute resolution
chol_indices = [1:4, 10];
chol_colors = lines(5); % Single color from MATLAB's "lines" colormap

% Define line styles for each treatment model
line_styles = {'-', '--', ':'}; % No drug, statin, ezetimibe

[t_PUFA, y_PUFA] = ode45(@(t, y) nodrugPUFAODE(t, y, p_PUFA), tspan, y0_PUFA);

% Run simulation for low-fat intake
[t, y] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0);

% Plot results
figure(1); clf; hold on;
models = {y, y_PUFA};
for i = 1:2
    for j = 1:length(chol_indices) % Iterate over compartments (GI, peripheral, liver, clearance)
        plot(t, models{i}(:, chol_indices(j)), ...
            'LineStyle', line_styles{i}, ...
            'Color', chol_colors(j, :), ...
            'LineWidth', 1.5);
    end
end

% plot(t_PUFA, y_PUFA(:, [1:4, 10]), 'LineWidth', 1.5); hold on;
% plot(t, y(:, [1:4, 10]), '--', 'LineWidth', 1.5);
title('Comparison of Saturated and Polyunsaturated Fatty Acid Intake on Lipid Dynamics');
xlabel('Time (hours)');
ylabel('Concentration (mg/dL)');
legend(...
       'SA Fat - GI Chol', 'SA Fat - Peripheral Chol', 'SA Fat - Liver Chol', ...
       'SA Fat - Clearance Chol', 'SA Fat - oxLDL', ...
       'PU Fat - GI Chol', 'PU Fat - Peripheral Chol', 'PU Fat - Liver Chol', ...
       'PU Fat - Clearance Chol', 'PU Fat - oxLDL' ...
       );
grid on;

