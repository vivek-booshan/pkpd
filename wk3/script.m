%% IV Compare Solution Methods
clf;

q = 0;
V = 1;
kcl = 1;
infusion = Model(q, V, kcl);
y0 = [1, 0];
tspan = 0:1/60:10;

[t, y] = ode45(@Model.infusionODE, tspan, y0, [], infusion.params);

figure(4); tiledlayout(1, 2)
% First plot: concentration-time profile
nexttile; hold on;
plot(t, exp(-t), 'b-', 'LineWidth', 2);  % Analytic solution
plot(t, y(:, 1), 'r--', 'LineWidth', 2);  % ODE45 solution
[teuler, yeuler] = Model.euler(@(t, y) Model.infusionODE(t, y, infusion.params), [0, 10], [1, 0]);
plot(teuler, yeuler(:, 1), 'g-.', 'LineWidth', 2);  % Euler solution

% Adding labels and legends
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
% title('IV Infusion: Concentration-Time Profile');
legend('Analytic Solution', 'ODE45 Solution', 'Euler Solution', 'Location', 'best');
grid on;

% Second plot: residual errors
nexttile; hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 1);  % Zero line
plot(t, y(:, 1) - exp(-t), 'r--', 'LineWidth', 2);  % Residuals ODE45
plot(teuler, yeuler(:, 1) - exp(-teuler), 'g-.', 'LineWidth', 2);  % Residuals Euler

% Adding labels and legends
xlabel('Time (hours)');
ylabel('Residuals (mg/L)');
title('Residuals: ODE vs Analytic Solution');
legend('Zero Line', 'ODE45 Residuals', 'Euler Residuals', 'Location', 'best');
grid on;

%% IV Dosing Plotting
figure(1);
% Define pastel RGB colors
colors = [0.5, 0.7, 0.9;  % Pastel blue
          0.7, 0.9, 0.5;  % Pastel green
          0.9, 0.7, 0.5;  % Pastel orange
          0.9, 0.5, 0.7]; % Pastel pink

% Define parameters
% q = 0; % mmol/hr
V = 1; % liters
kcl = 1; % 1/hr
ka = 1; % 1/hr
oral = Model(0, V, kcl);
y0 = [1, 0];
tspan = 0:1/60:10; % Define time span for integration

% Noise functions
error_prop = @(samples) 0.05 * samples .* randn(size(samples));
error_add = @(samples) 0.05 * mean(samples) .* randn(size(samples));
analytic = @(params, t) params(1)*exp(-params(2)*t);
lsq_options = optimset('Display', 'off');
param_fit = @(t, y) lsqcurvefit(analytic, randn(1, 2), t, y, [], [], lsq_options);

% Solve ODE
[t, y] = ode45(@Model.infusionODE, tspan, y0, [], oral.params);

% Create figure with tiled layout
figure('Color', [0.95, 0.95, 0.9]); % Cream background color
tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% Loop through scenarios
for i = 0:3
    % Setup for each tile
    nexttile; hold on;
    %set(gca, 'Color', [0.95, 0.95, 0.9])
    % Define slice and samples for the scenario
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);
    
    % Basic Fit
    param_est_basic = param_fit(t_samples, samples);
    fitted_eqn_basic = @(t) analytic(param_est_basic, t);
    plot(t_samples, samples, 'o', 'MarkerFaceColor', colors(1, :), 'DisplayName', 'Basic Samples');
    plot(tspan, fitted_eqn_basic(tspan), '-', 'LineWidth', 2, 'Color', colors(1, :), 'DisplayName', 'Basic Fitted Curve');
    
    % Proportional Noise Fit
    noisy_prop = samples + error_prop(samples);
    param_est_prop = param_fit(t_samples, noisy_prop);
    fitted_eqn_prop = @(t) analytic(param_est_prop, t);
    plot(t_samples, noisy_prop, 'o', 'MarkerFaceColor', colors(3, :), 'DisplayName', 'Proportional Samples');
    plot(tspan, fitted_eqn_prop(tspan), '--', 'LineWidth', 2, 'Color', colors(3, :), 'DisplayName', 'Proportional Fitted Curve');
    
    % Additive Noise Fit
    noisy_add = samples + error_add(samples);
    param_est_add = param_fit(t_samples, noisy_add);
    fitted_eqn_add = @(t) analytic(param_est_add, t);
    plot(t_samples, noisy_add, 'o', 'MarkerFaceColor', colors(2, :), 'DisplayName', 'Additive Samples');
    plot(tspan, fitted_eqn_add(tspan), '-.', 'LineWidth', 2, 'Color', colors(2, :), 'DisplayName', 'Additive Fitted Curve');
    
    % Calculate R^2 for each fit
    R_squared_basic = Model.R2(samples, fitted_eqn_basic(t_samples));
    R_squared_prop = Model.R2(noisy_prop, fitted_eqn_prop(t_samples));
    R_squared_add = Model.R2(noisy_add, fitted_eqn_add(t_samples));
    
    % Add R^2 annotations
    text(10, 0.9, sprintf('Basic R^2 = %f', R_squared_basic), 'FontSize', 12, 'FontWeight', 'bold', 'Color', colors(4, :), 'HorizontalAlignment', 'right');
    text(10, 0.8, sprintf('Proportional R^2 = %f', R_squared_prop), 'FontSize', 12, 'FontWeight', 'bold', 'Color', colors(4, :), 'HorizontalAlignment', 'right');
    text(10, 0.7, sprintf('Additive R^2 = %f', R_squared_add), 'FontSize', 12, 'FontWeight', 'bold', 'Color', colors(4, :), 'HorizontalAlignment', 'right');
    
    % Label axes and add title
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Concentration (mmol/L)', 'FontSize', 12);
    title(sprintf('%d Samples', 4-i), 'FontSize', 14, 'FontWeight', 'bold');
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    
    % Add legend
    legend('show', 'Location', 'north', 'FontSize', 10);
end

%% Oral Dosing Plotting
figure(2);
% Define parameters
q = 0; % mmol/hr
V = 1; % liters
kcl = 0.01; % 1/hr
ka = 0.4; % 1/hr
oral = Model(q, V, kcl, ka);
y0 = [1, 0, 0];
tspan = 0:1/60:10; % Define time span for integration

% Noise functions
error_prop = @(samples) 0.05 * samples .* randn(size(samples));
error_add = @(samples) 0.05 * mean(samples) .* randn(size(samples));
analytic = @(params, t) (params(1)/(params(2) - params(1))) * (exp(-params(1)*t) - exp(-params(2)*t));
lsq_options = optimset('Display', 'off');
param_fit = @(t, y) lsqcurvefit(analytic, [1, 0, 0], t, y, [], [], lsq_options);

% Solve ODE
[t, y] = ode45(@Model.oralODE, tspan, y0, [], oral.params);

% Define colors
colors = [0.5, 0.7, 0.9;  % Pastel blue
          0.7, 0.9, 0.5;  % Pastel green
          0.9, 0.7, 0.5;  % Pastel orange
          0.9, 0.5, 0.7]; % Pastel pink

% Create figure with tiled layout
figure('Color', [0.95, 0.95, 0.9]); % Cream background for the figure
tiledlayout(2, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Loop through scenarios
for i = 0:3
    % Setup for each tile
    nexttile; hold on;
    %set(gca, 'Color', [0.95, 0.95, 0.9]); % Cream background for the tiles
    
    % Define slice and samples for the scenario
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);
    
    % Basic Fit
    param_est_basic = param_fit(t_samples, samples);
    fitted_eqn_basic = @(t) analytic(param_est_basic, t);
    plot(t_samples, samples, 'o', 'MarkerFaceColor', colors(1, :), 'MarkerEdgeColor', colors(1, :), 'DisplayName', 'Basic Samples');
    plot(tspan, fitted_eqn_basic(tspan), '-', 'LineWidth', 2, 'Color', colors(1, :), 'DisplayName', 'Basic Fitted Curve');
    
    % Proportional Noise Fit
    noisy_prop = samples + error_prop(samples);
    param_est_prop = param_fit(t_samples, noisy_prop);
    fitted_eqn_prop = @(t) analytic(param_est_prop, t);
    plot(t_samples, noisy_prop, 's', 'MarkerFaceColor', colors(2, :), 'MarkerEdgeColor', colors(2, :), 'DisplayName', 'Proportional Samples');
    plot(tspan, fitted_eqn_prop(tspan), '--', 'LineWidth', 2, 'Color', colors(2, :), 'DisplayName', 'Proportional Fitted Curve');
    
    % Additive Noise Fit
    noisy_add = samples + error_add(samples);
    param_est_add = param_fit(t_samples, noisy_add);
    fitted_eqn_add = @(t) analytic(param_est_add, t);
    plot(t_samples, noisy_add, 'd', 'MarkerFaceColor', colors(3, :), 'MarkerEdgeColor', colors(1, :), 'DisplayName', 'Additive Samples');
    plot(tspan, fitted_eqn_add(tspan), '-.', 'LineWidth', 2, 'Color', colors(3, :), 'DisplayName', 'Additive Fitted Curve');
    
    % Calculate R^2 for each fit
    R_squared_basic = Model.R2(samples, fitted_eqn_basic(t_samples));
    R_squared_prop = Model.R2(noisy_prop, fitted_eqn_prop(t_samples));
    R_squared_add = Model.R2(noisy_add, fitted_eqn_add(t_samples));
    
    % Add R^2 annotations
    text(10, 0.9, sprintf('Basic R^2 = %f', R_squared_basic), 'FontSize', 12, 'FontWeight', 'bold', 'Color', colors(4, :), 'HorizontalAlignment', 'right');
    text(10, 0.8, sprintf('Proportional R^2 = %f', R_squared_prop), 'FontSize', 12, 'FontWeight', 'bold', 'Color', colors(4, :), 'HorizontalAlignment', 'right');
    text(10, 0.7, sprintf('Additive R^2 = %f', R_squared_add), 'FontSize', 12, 'FontWeight', 'bold', 'Color', colors(4, :), 'HorizontalAlignment', 'right');
    
    % Label axes and add title
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Concentration (mmol/L)', 'FontSize', 12);
    title(sprintf('%d Samples', 4-i), 'FontSize', 14, 'FontWeight', 'bold');
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    
    % Add legend
    legend('show', 'Location', 'north', 'FontSize', 10);
end

%% Compare IV and Oral Dosing
dose = [1, 1, 1];
V = [1, 1, 1];
ka = [1, 1, 0.5];
kcl = [1, 0.5, 1];
tspan = 0:1/60:10; % Define time span for integration

% Define a modern color palette similar to Seaborn's HSL style
colors = [0.1, 0.2, 0.5; % Deep blue
          0.9, 0.3, 0.3; % Bright red
          0.3, 0.7, 0.5]; % Soft green

figure('Color', 'w'); % Set figure background color
tiledlayout(1, 3, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for i = 1:3
    % Setup infusion and oral models
    infusion = Model(0, V(i), kcl(i));
    oral = Model(0, V(i), kcl(i), ka(i));
    
    % Initial conditions
    y0_infusion = [1, 0];
    y0_oral = [dose(i), 0, 0];
    
    % Solve ODEs
    [t, y1] = ode45(@Model.infusionODE, tspan, y0_infusion, [], infusion.params);
    [~, y2] = ode45(@Model.oralODE, tspan, y0_oral, [], oral.params);
    
    % Plotting
    nexttile; hold on;
    
    plot(t, y1(:, 1), 'Color', colors(1, :), 'LineWidth', 2, 'DisplayName', 'IV', 'MarkerSize', 8);
    plot(t, y2(:, 2), 'Color', colors(2, :), 'LineWidth', 2, 'DisplayName', 'Oral', 'MarkerSize', 8);
    
    % Add labels and titles
    xlabel('Time (hours)', 'FontSize', 12);
    ylabel('Concentration (mmol/L)', 'FontSize', 12);
    title(sprintf('Condition %d', i), 'FontSize', 14, 'FontWeight', 'bold');
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    ax = gca;
    ax.XColor = [0.2, 0.2, 0.2];
    ax.YColor = [0.2, 0.2, 0.2];
    
    % Adjust legend position
    legend('show', 'Location', 'northeast', 'FontSize', 10);
end

%% Oral compare integrators
figure(5);
q = 0; % mmol/hr
V = 1; % liters
kcl = 0.5; % 1/hr
ka = 1; % 1/hr
oral = Model(q, V, kcl, ka);
y0 = [1, 0, 0];
tspan = 0:1/60:10; % Define time span for integration

% Noise functions
error_prop = @(samples) 0.05 * samples .* randn(size(samples));
error_add = @(samples) 0.05 * mean(samples) .* randn(size(samples));
analytic = @(t) (ka/(kcl - ka))*(exp(-ka*t) - exp(-kcl*t));
% Solve ODE
[t, y] = ode45(@Model.oralODE, tspan, y0, [], oral.params);

figure(4); tiledlayout(1, 2)
% First plot: concentration-time profile
nexttile; hold on;
plot(t, analytic(t), 'b-', 'LineWidth', 2);  % Analytic solution
plot(t, y(:, 2), 'LineWidth', 2);  % ODE45 solution
[teuler, yeuler] = Model.euler(@(t, y) Model.oralODE(t, y, oral.params), [0, 10], [1, 0, 0]);
plot(teuler, yeuler(:, 2), 'g-.', 'LineWidth', 2);  % Euler solution

% Adding labels and legends
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('IV Infusion: Concentration-Time Profile');
legend('Analytic Solution', 'ODE45 Solution', 'Euler Solution', 'Location', 'best');
grid on;

% Second plot: residual errors
nexttile; hold on;
plot(t, zeros(size(t)), 'k--', 'LineWidth', 1);  % Zero line
plot(t, y(:, 2) - analytic(t), 'r--', 'LineWidth', 2);  % Residuals ODE45
plot(teuler, yeuler(:, 2) - analytic(teuler), 'g-.', 'LineWidth', 2);  % Residuals Euler

% Adding labels and legends
xlabel('Time (hours)');
ylabel('Residuals (mg/L)');
title('Residuals: ODE vs Analytic Solution');
legend('Zero Line', 'ODE45 Residuals', 'Euler Residuals', 'Location', 'best');
grid on;
