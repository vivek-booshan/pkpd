params.KTR_LCIG = 9.2;       % (1/hr) for LCIG
params.KTR_LC_oral = 2.4;    % (1/hr) for LC-oral
params.CL_F = 24.8;          % Clearance (L/hr)
params.Vc_F = 58.5;          % Central volume (L)
params.Q_F = 6.8;            % Inter-compartmental clearance (L/hr)
params.Vp_F = 72.9;          % Peripheral volume (L)
params.F_rel_LCIG = 0.97;    % Bioavailability for LCIG relative to LC-oral
params.formulation = 'LCIG'; % Use 'LCIG' or 'LC-oral'

% Initial conditions
dose = 100;
y0 = [dose; 0; 0; 0];  % Assume dose in the gut and 0 in other compartments initially

% Time span (e.g., 0 to 24 hours)
tspan = [0 24];

% Solve ODE
[t, y] = ode45(@(t, y) levodopa_ODE(t, y, params), tspan, y0);
balance = 100 - sum(y, 2);  % Calculate balance

% Create a tiled layout for the plots
figure;
tiledlayout(2, 1, 'Padding', 'compact'); % Compact layout

% Plot drug amounts in each compartment
nexttile;
hold on; % Allow multiple plots on the same axis
plot(t, y(:,1), 'LineWidth', 2, 'Color', [0.2, 0.4, 0.8]);  % Gut
plot(t, y(:,2), 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2]);  % Central
plot(t, y(:,3), 'LineWidth', 2, 'Color', [0.2, 0.8, 0.2]);  % Peripheral
plot(t, y(:,4), 'LineWidth', 2, 'Color', [0.6, 0.4, 0.8]);  % Lost
hold off;

% Add titles, labels, and legend
title('Drug Amounts in Compartments', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Amount (mg)', 'FontSize', 12);
legend({'Gut', 'Central', 'Peripheral', 'Lost'}, 'Location', 'northeast', 'FontSize', 10);
grid on; % Add grid for better readability
set(gca, 'FontSize', 12); % Set axis font size

% Plot balance over time
nexttile;
plot(t, balance, 'LineWidth', 2, 'Color', [0.9, 0.6, 0], 'LineStyle', '--'); % Use a dashed line for balance
title('Balance of Drug Amount', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Balance (mg)', 'FontSize', 12);
grid on; % Add grid for better readability
set(gca, 'FontSize', 12); % Set axis font size

% Overall figure adjustments
sgtitle('Pharmacokinetic Model Results', 'FontSize', 16, 'FontWeight', 'bold');
