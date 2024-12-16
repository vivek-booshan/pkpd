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
    'clearance', struct('chol', 0.01, 'TG', 0.01) ...
);

% Liver parameters (transport and clearance rates, assumed in 1/hr)
p.liver = struct( ...
    'LDL', struct('chol', 0.07, 'TG', 0.05), ...
    'HDL', struct('chol', 0.04, 'TG', 0.03), ...
    'clearance', struct('chol', 0.02, 'TG', 0.01), ...
    'oxLDL', struct('chol', 0.005) ...
);

% ROS dynamics parameters
p.basalROS = 0.02;                % Baseline ROS production (arbitrary units/hr)
p.foodProductionMultiplier = 1.2; % ROS increase due to food intake (unitless multiplier)
p.antioxidant = 0.01;             % Antioxidant clearance of ROS (1/hr)
p.foodAntiOxMultiplier = 1.0;     % Antioxidant effect (unitless multiplier)

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
    0       % oxLDL (starts at 0)
]';
%%

tspan = 1:1/60:24;
[t, y] = ode45(@(t, y) nodrugODE(t, y, p), tspan, y0);
% volumes = [p.V.GI, p.V.peripheral, p.V.liver, 1];
% b = sum(volumes)
figure(1); clf
plot(y(:, [1:4 10]))
legend()
