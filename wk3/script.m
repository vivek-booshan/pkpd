%% IV Compare Solution Methods
clf;

q = 0;
V = 1;
kcl = 1;
infusion = Model(q, V, kcl);
y0 = [1, 0];
tspan = 0:1/60:10;

error_prop = @(samples) 0.05 * samples .* randn(size(samples));
error_add = @(samples) 0.05 * mean(samples) .* randn(size(samples));
analytic = @(params, t) params(1)*exp(-params(2)*t);
lsq_options = optimset('Display', 'off');
param_fit = @(t, y) lsqcurvefit(analytic, [randn(1), randn(1)], t, y, [], [], lsq_options);
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
title('IV Infusion: Concentration-Time Profile');
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
%% IV fit
% Basic fit
figure(1);
tiledlayout(2, 2);
disp("Basic fit")
for i = 0:3
    nexttile; hold on;
    
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');
    
    param_est = param_fit(t_samples, samples);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    R_squared = Model.R2(samples, fitted_eqn(t_samples));
    
    % Add R^2 annotation
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');

    % Label axes
    xlabel('Time (hours)');
    ylabel('Concentration (mg/L)');
    title(sprintf('Basic Fit - Scenario %d', i+1));
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    legend('show', 'Location', 'east');
end

% Proportional noise fit
figure(2);
disp("\n Proportional fit")
tiledlayout(2, 2);
for i = 0:3
    nexttile; hold on;
    
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);
    noisy_prop = samples + error_prop(samples);

    plot(t_samples, noisy_prop, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');
    
    param_est = param_fit(t_samples, noisy_prop);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    R_squared = Model.R2(noisy_prop, fitted_eqn(t_samples));
    
    % Add R^2 annotation
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');

    % Label axes
    xlabel('Time (hours)');
    ylabel('Concentration (mg/L)');
    title(sprintf('Proportional Noise Fit - Scenario %d', i+1));
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    legend('show', 'Location', 'east');
end

% Additive noise fit
figure(3);
tiledlayout(2, 2);
disp("\n Additive fit");
for i = 0:3
    nexttile; hold on;
    
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    
    noisy_add = samples + error_add(samples);
    plot(t_samples, samples, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');

    param_est = param_fit(t_samples, noisy_add);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    R_squared = Model.R2(noisy_add, fitted_eqn(t_samples));
    
    % Add R^2 annotation
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');

    % Label axes
    xlabel('Time (hours)');
    ylabel('Concentration (mg/L)');
    title(sprintf('Additive Noise Fit - Scenario %d', i+1));
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    legend('show', 'Location', 'east');
end

%% Oral Dosing Plotting

q = 0; % mmol/hr
V = 1; % liters
kcl = 0.01; % 1/hr
ka = 0.4; % 1/hr
oral = Model(q, V, kcl, ka);
y0 = [1, 0, 0];
tspan = 0:1/60:10; % Define time span for integration
error_prop = @(samples) 0.05 * samples .* randn(size(samples));
error_add = @(samples) 0.05 * mean(samples) .* randn(size(samples));
analytic = @(params, t) (params(1)/(params(2) + params(1)))*(exp(-params(2)*t) - exp(-params(1)*t));
lsq_options = optimset('Display', 'off');
param_fit = @(t, y) lsqcurvefit(analytic, [1, 0, 0], t, y, [], [], lsq_options);

% Solve ODE
[t, y] = ode45(@Model.oralODE, tspan, y0, [], oral.params);

% Basic fit
figure(1);
tiledlayout(2, 2);
disp("Basic fit");
for i = 0:3
    nexttile; hold on;
    
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');
    
    param_est = param_fit(t_samples, samples);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    R_squared = Model.R2(samples, fitted_eqn(t_samples));
    
    % Add R^2 annotation
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');

    % Label axes
    xlabel('Time (hours)');
    ylabel('Concentration (mmol/L)');
    title(sprintf('Basic Fit - Scenario %d', i+1));
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;
    
    % Move legend outside the plot
    legend('show', 'Location', 'east');
end

% Fit proportional noise sample
figure(2);
disp("\n Proportional fit");
tiledlayout(2, 2);
for i = 0:3
    nexttile; hold on;
    
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');
    
    noisy_prop = samples + error_prop(samples);
    param_est = param_fit(t_samples, noisy_prop);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    R_squared = Model.R2(noisy_prop, fitted_eqn(t_samples));
    
    % Add R^2 annotation
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');

    % Label axes
    xlabel('Time (hours)');
    ylabel('Concentration (mmol/L)');
    title(sprintf('Proportional Noise Fit - Scenario %d', i+1));
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;

    % Move legend outside the plot
    legend('show', 'Location', 'east');
end

% Fit additive noise sample
figure(3);
tiledlayout(2, 2);
disp("\n Additive fit");
for i = 0:3
    nexttile; hold on;
    
    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Samples');
    
    noisy_add = samples + error_add(samples);
    param_est = param_fit(t_samples, noisy_add);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'b-', 'LineWidth', 2, 'DisplayName', 'Fitted Curve');
    R_squared = Model.R2(noisy_add, fitted_eqn(t_samples));
    
    % Add R^2 annotation
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');

    % Label axes
    xlabel('Time (hours)');
    ylabel('Concentration (mmol/L)');
    title(sprintf('Additive Noise Fit - Scenario %d', i+1));
    xlim([0, 10]);
    ylim([0, 1]);
    grid on;

    % Move legend outside the plot
    legend('show', 'Location', 'east');
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
