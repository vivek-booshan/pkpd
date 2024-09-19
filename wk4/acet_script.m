clf;
% Define parameters
p = acet_param();
p.KEP = 0.07;
p.KMG = 0.1;
p.VMG = 50;
p.KEG = 0.2;
p.KMS = 0.05;
p.VMS = 40;
p.KES = 0.3;
p.VMCM = 10;
p.FEP = 0.1;
p.KMCM = 0.05;
p.DOSE1 = 60;  % Example dose in mg/kg
p.DOSE2 = 90;  % Example dose in mg/kg
p.TLAG1 = p.tlag;    % Time lag for dose1 in hours
p.TLAG2 = p.tlag;    % Time lag for dose2 in hours
p.V = 70;       % Volume of distribution in L
p.KA = 1;       % Absorption rate constant in 1/hr

% Initial conditions
y0 = zeros(14, 1);

% Time span for simulation
tspan = [0 24];  % Simulate for 24 hours

% Solve the ODE
[t, y] = ode45(@(t, y) acet_ODE(t, y, p), tspan, y0);

% Define some variables to beautify the plot
alpha = 0.5; % Alpha transparency for second dose
colors = lines(4); % Use the lines color palette for distinct colors
thicknessFirstDose = 2; % Line thickness for the first dose
thicknessSecondDose = 1; % Line thickness for the second dose

% Create a figure with better layout
figure(1);

% Plasma Concentrations
subplot(2, 1, 1);
hold on;
plot(t, y(:, 1), '-', 'DisplayName', 'Cp (1st Dose)', 'Color', colors(1, :), 'LineWidth', thicknessFirstDose);
plot(t, y(:, 2), '-', 'DisplayName', 'CpG (1st Dose)', 'Color', colors(2, :), 'LineWidth', thicknessFirstDose);
plot(t, y(:, 3), '-', 'DisplayName', 'CpS (1st Dose)', 'Color', colors(3, :), 'LineWidth', thicknessFirstDose);

plot(t, y(:, 4), '--', 'DisplayName', 'Cp (2nd Dose)', 'Color', [colors(1, :) alpha], 'LineWidth', thicknessSecondDose);
plot(t, y(:, 5), '--', 'DisplayName', 'CpG (2nd Dose)', 'Color', [colors(2, :) alpha], 'LineWidth', thicknessSecondDose);
plot(t, y(:, 6), '--', 'DisplayName', 'CpS (2nd Dose)', 'Color', [colors(3, :) alpha], 'LineWidth', thicknessSecondDose);

xlabel('Time (hours)', 'FontSize', 12);
ylabel('Concentration (mg/L)', 'FontSize', 12);
title('Plasma Concentrations of Acetaminophen and Metabolites', 'FontSize', 14);
legend('Location', 'best');
grid on; % Optional: Add grid for better visibility

% Cumulative Urine Recoveries
subplot(2, 1, 2);
hold on;
plot(t, y(:, 7), '-', 'DisplayName', 'UP (1st Dose)', 'Color', colors(1, :), 'LineWidth', thicknessFirstDose);
plot(t, y(:, 8), '-', 'DisplayName', 'UG (1st Dose)', 'Color', colors(2, :), 'LineWidth', thicknessFirstDose);
plot(t, y(:, 9), '-', 'DisplayName', 'US (1st Dose)', 'Color', colors(3, :), 'LineWidth', thicknessFirstDose);
plot(t, y(:, 10), '-', 'DisplayName', 'UCM (1st Dose)', 'Color', colors(4, :), 'LineWidth', thicknessFirstDose);

plot(t, y(:, 11), '--', 'DisplayName', 'UP (2nd Dose)', 'Color', [colors(1, :) alpha], 'LineWidth', thicknessSecondDose);
plot(t, y(:, 12), '--', 'DisplayName', 'UG (2nd Dose)', 'Color', [colors(2, :) alpha], 'LineWidth', thicknessSecondDose);
plot(t, y(:, 13), '--', 'DisplayName', 'US (2nd Dose)', 'Color', [colors(3, :) alpha], 'LineWidth', thicknessSecondDose);
plot(t, y(:, 14), '--', 'DisplayName', 'UCM (2nd Dose)', 'Color', [colors(4, :) alpha], 'LineWidth', thicknessSecondDose);

xlabel('Time (hours)', 'FontSize', 12);
ylabel('Cumulative Recovery (mg)', 'FontSize', 12);
title('Cumulative Urine Recoveries of Acetaminophen Metabolites', 'FontSize', 14);
legend('Location', 'best');
grid on; % Optional: Add grid for better visibility

% Set figure properties for better clarity
set(gcf, 'Color', 'w');