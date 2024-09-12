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
%% IV Fits
% basic fit
figure(1);
tiledlayout(2, 2)
disp("Basic fit")
for i = 0:3
    nexttile; hold on;

    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro');
    
    param_est = param_fit(t_samples, samples);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'r');
    R_squared = Model.R2(samples, fitted_eqn(t_samples));
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');
    xlim([0, 10]);
    ylim([0, 1]);
end

% fit proportional noise sample
figure(2);
disp("\n Proportional fit")
tiledlayout(2, 2);
for i = 0:3
    nexttile; hold on;

    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro');
    
    noisy_prop = samples + error_prop(samples);
    % lsq_options = optimset('Display', 'off');
    param_est = param_fit(t_samples, noisy_prop);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'r');
    R_squared = Model.R2(noisy_prop, fitted_eqn(t_samples));
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');
    xlim([0, 10]);
    ylim([0, 1]);
    %disp(Model.R2(samples, fitted_eqn(t_samples)));
end

% fit additive noise sample
figure(3);
tiledlayout(2, 2);
disp("\n Additive fit");
for i = 0:3
    nexttile; hold on;

    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro');
    
    noisy_add = samples + error_add(samples);

    % lsq_options = optimset('Display', 'off');
    param_est = param_fit(t_samples, noisy_add);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'r');
    R_squared = Model.R2(noisy_add, fitted_eqn(t_samples));
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');
    xlim([0, 10]);
    ylim([0, 1]);
end

%% Oral
q = 0; % mmol/hr
V = 1; % liters
kcl = 0.04; % 1/hr
ka = 0.2; % mmol/hr
oral = Model(q, V, kcl, ka);
y0 = [1, 0, 0];
%y0 = [oral.params.q, 0, 0];
error_prop = @(samples) 0.05 * samples .* randn(size(samples));
error_add = @(samples) 0.05 * mean(samples) .* randn(size(samples));
analytic = @(params, t) (params(1)/(params(2) + params(1)))*(exp(-params(2)*t) - exp(-params(1)*t));
lsq_options = optimset('Display', 'off');
param_fit = @(t, y) lsqcurvefit(analytic, [1, 0, 0], t, y, [], [], lsq_options);

[t, y] = ode45(@Model.oralODE, tspan, y0, [], oral.params);
figure(1);
tiledlayout(2, 2)
disp("Basic fit")
for i = 0:3
    nexttile; hold on;

    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro');
    
    lsq_options = optimset('Display', 'off');
    param_est = param_fit(t_samples, samples);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'r');
    R_squared = Model.R2(samples, fitted_eqn(t_samples));
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');
    xlim([0, 10]);
    ylim([0, 1]);
end

% fit proportional noise sample
figure(2);
disp("\n Proportional fit")
tiledlayout(2, 2);
for i = 0:3
    nexttile; hold on;

    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro');
    
    noisy_prop = samples + error_prop(samples);
    lsq_options = optimset('Display', 'off');
    param_est = param_fit(t_samples, noisy_prop);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'r');
    R_squared = Model.R2(noisy_prop, fitted_eqn(t_samples));
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');
    xlim([0, 10]);
    ylim([0, 1]);
    %disp(Model.R2(samples, fitted_eqn(t_samples)));
end

% fit additive noise sample
figure(3);
tiledlayout(2, 2);
disp("\n Additive fit");
for i = 0:3
    nexttile; hold on;

    slice = 30 + 120*(0:4-i);
    samples = y(slice, 1);
    t_samples = t(slice);

    plot(t_samples, samples, 'ro');
    
    noisy_add = samples + error_add(samples);

    lsq_options = optimset('Display', 'off');
    param_est = param_fit(t_samples, noisy_add);
    fitted_eqn = @(t) analytic(param_est, t);

    plot(tspan, fitted_eqn(tspan), 'r');
    R_squared = Model.R2(noisy_add, fitted_eqn(t_samples));
    text(9, 0.9, sprintf('R^2 = %f', R_squared), 'FontSize', 12, 'HorizontalAlignment', 'right');
    xlim([0, 10]);
    ylim([0, 1]);
end


%% compare
%q = [1, 1, 1];
dose = [1, 1, 1];
V = [1, 1, 1];
ka = [1, 1, 0.5];
kcl = [1, 0.5, 1];

tiledlayout(1, 3);
for i = 1:3
    infusion = Model(0, V(i), kcl(i));
    oral = Model(0, V(i), kcl(i), ka(i));
    y0_infusion = [1, 0];
    y0_oral = [dose(i), 0, 0];
    [~, y1] = ode45(@Model.infusionODE, tspan, y0_infusion, [], infusion.params);
    [t, y2] = ode45(@Model.oralODE, tspan, y0_oral, [], oral.params);
    
    nexttile; hold on;
    plot(t, y1(:, 1)); % (hr, mmol/hr)
    plot(t, y2(:, 2)); % (hr, mmol/hr)
    legend('IV', 'oral')
end
