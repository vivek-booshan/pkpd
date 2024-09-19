% Parameters
k_e = 0.1; % Elimination rate constant (1/h)
V_d = 10; % Volume of distribution (L)
U = 0.5; % Concentration of drug in urine (mg/L)
V = 1; % Urine flow rate (L/h)
C_p_initial = 10; % Initial plasma concentration (mg/L)
Q = 5; % Liver blood flow (L/h)
C_in = 5; % Concentration entering the liver (mg/L)
C_out = 1; % Concentration leaving the liver (mg/L)
V_max = 0.8; % Maximum rate of metabolism (mg/h)
K_m = 0.5; % Michaelis-Menten constant (mg/L)

% Time vector
t = 0:0.1:24; % Time from 0 to 24 hours

% First Order Clearance Profile
C_first_order = C_p_initial * exp(-k_e * t); % Concentration over time with first-order elimination

% Renal Clearance Profile
% Assuming renal clearance directly correlates to first-order elimination for simplicity
C_renal = C_p_initial * exp(-k_e * t); % Same as first-order for simplicity

% Hepatic Clearance Profile
% To reflect hepatic processing, assuming an exponential decrease
C_hepatic = C_p_initial * exp(-Q / V_d * t); % Concentration over time with hepatic clearance

% Michaelis-Menten Clearance Profile
% Solve the differential equation using ode45
[t_michaelis, C_michaelis] = ode45(@(t, Cp) -V_max * Cp / (K_m + Cp), t, C_p_initial);

% Plotting
figure;
hold on;

% Define HSL colors
colors = [hsl2rgb(180, 0.6, 0.6); % First Order Clearance (light blue)
          hsl2rgb(120, 0.6, 0.6); % Renal Clearance (light green)
          hsl2rgb(0, 0.6, 0.6);   % Hepatic Clearance (light red)
          hsl2rgb(240, 0.6, 0.6)]; % Michaelis-Menten Clearance (light purple)

% Plot each profile with the custom colors
%plot(t, C_first_order, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'First Order Clearance');
plot(t, C_renal, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'Renal Clearance');
plot(t, C_hepatic, 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', 'Hepatic Clearance');
%plot(t_michaelis, C_michaelis, 'Color', colors(4,:), 'LineWidth', 2, 'DisplayName', 'Michaelis-Menten Clearance');

hold off;

% Formatting the plot
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
xlabel('Time (hours)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Concentration (mg/L)', 'FontSize', 14, 'FontWeight', 'bold');
title('Drug Concentration Depletion with Different Clearance Mechanisms', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'northeast');
grid on;
xlim([0 24]);
ylim([0 max([C_first_order, C_renal, C_hepatic, C_michaelis]) * 1.1]); % Adjust y-axis to fit all profiles

% Function to convert HSL to RGB
function rgb = hsl2rgb(h, s, l)
    c = (1 - abs(2 * l - 1)) * s;
    x = c * (1 - abs(mod(h / 60, 2) - 1));
    m = l - c / 2;
    switch floor(h / 60)
        case 0
            rgb = [c, x, 0] + m;
        case 1
            rgb = [x, c, 0] + m;
        case 2
            rgb = [0, c, x] + m;
        case 3
            rgb = [0, x, c] + m;
        case 4
            rgb = [x, 0, c] + m;
        case 5
            rgb = [c, 0, x] + m;
    end
end

% Clear workspace
clear;
clc;

% Parameters
params = struct();

% Initial plasma concentration
C_p_initial = 10; % Initial plasma concentration (mg/L)
t = 0:0.1:24; % Time vector from 0 to 24 hours

% Define parameters for each clearance mechanism

% Renal Clearance Parameters
params.renal.k_e = 0.1; % Elimination rate constant (1/h)
params.renal.U = 0.5; % Concentration of drug in urine (mg/L)
params.renal.V = 1; % Urine flow rate (L/h)

% Hepatic Clearance Parameters
params.hepatic.Q = 5; % Liver blood flow (L/h)
params.hepatic.V_d = 10; % Volume of distribution (L)

% Michaelis-Menten Clearance Parameters
params.michaelis.V_max = 0.8; % Maximum rate of metabolism (mg/h)
params.michaelis.K_m = 0.5; % Michaelis-Menten constant (mg/L)

% Define ranges for parameter variations
renal_k_e_range = [0.05, 0.1, 0.2]; % Variation in elimination rate constant (1/h)
hepatic_Q_range = [2, 5, 10]; % Variation in liver blood flow (L/h)
hepatic_V_d_range = [5, 10, 20]; % Variation in volume of distribution (L)
michaelis_V_max_range = [0.4, 0.8, 1.2]; % Variation in maximum rate of metabolism (mg/h)
michaelis_K_m_range = [0.3, 0.5, 0.7]; % Variation in Michaelis-Menten constant (mg/L)

% Create figure for Renal Clearance
figure;
subplot(1,1,1); % Single subplot
hold on;
for k_e = renal_k_e_range
    C_renal = C_p_initial * exp(-k_e * t);
    plot(t, C_renal, 'DisplayName', sprintf('k_e = %.2f 1/h', k_e));
end
hold off;
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('Renal Clearance');
legend('show');
grid on;

% Create figure for Hepatic and Michaelis-Menten Clearance
figure;

% Plot for Hepatic Clearance (varying liver blood flow)
subplot(2,2,1); % 2 rows, 2 columns, 1st subplot
hold on;
for Q = hepatic_Q_range
    C_hepatic = C_p_initial * exp(-Q / params.hepatic.V_d * t);
    plot(t, C_hepatic, 'DisplayName', sprintf('Q = %.1f L/h', Q));
end
hold off;
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('Hepatic Clearance - Varying Liver Blood Flow');
legend('show');
grid on;

% Plot for Hepatic Clearance (varying volume of distribution)
subplot(2,2,2); % 2 rows, 2 columns, 2nd subplot
hold on;
for V_d = hepatic_V_d_range
    C_hepatic_Vd = C_p_initial * exp(-params.hepatic.Q / V_d * t);
    plot(t, C_hepatic_Vd, 'DisplayName', sprintf('V_d = %.1f L', V_d));
end
hold off;
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('Hepatic Clearance - Varying Volume of Distribution');
legend('show');
grid on;

% Plot for Michaelis-Menten Clearance (varying V_max)
subplot(2,2,3); % 2 rows, 2 columns, 3rd subplot
hold on;
for V_max = michaelis_V_max_range
    [t_michaelis, C_michaelis] = ode45(@(t, Cp) -V_max * Cp / (params.michaelis.K_m + Cp), t, C_p_initial);
    plot(t_michaelis, C_michaelis, 'DisplayName', sprintf('V_{max} = %.1f mg/h', V_max));
end
hold off;
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('Michaelis-Menten Clearance - Varying V_{max}');
legend('show');
grid on;

% Plot for Michaelis-Menten Clearance (varying K_m)
subplot(2,2,4); % 2 rows, 2 columns, 4th subplot
hold on;
for K_m = michaelis_K_m_range
    [t_michaelis, C_michaelis] = ode45(@(t, Cp) -params.michaelis.V_max * Cp / (K_m + Cp), t, C_p_initial);
    plot(t_michaelis, C_michaelis, 'DisplayName', sprintf('K_m = %.1f mg/L', K_m));
end
hold off;
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('Michaelis-Menten Clearance - Varying K_m');
legend('show');
grid on;

% Adjust figure layout
sgtitle('Effects of Parameter Variations on Drug Clearance'); % Overall title for all subplots

% Clear workspace
clear;
clc;

% Parameters
k_e = 0.1; % Elimination rate constant (1/h)
V_d = 10; % Volume of distribution (L)
U = 0.5; % Concentration of drug in urine (mg/L)
V = 1; % Urine flow rate (L/h)
C_p_initial = 10; % Initial plasma concentration (mg/L)
Q = 5; % Liver blood flow (L/h)
C_in = 5; % Concentration entering the liver (mg/L)
C_out = 1; % Concentration leaving the liver (mg/L)
V_max = 0.8; % Maximum rate of metabolism (mg/h)
K_m = 0.5; % Michaelis-Menten constant (mg/L)
protein_binding = 0.9; % Fraction of drug bound to plasma proteins (0 to 1)

% Calculate free fraction of the drug
free_fraction = 1 - protein_binding;

% Adjust parameters for free drug
k_e_free = k_e * free_fraction; % Adjusted elimination rate constant
Q_free = Q * free_fraction; % Adjusted liver blood flow for free drug

% Time vector
t = 0:0.1:24; % Time from 0 to 24 hours

% First Order Clearance Profile
C_first_order = C_p_initial * exp(-k_e_free * t); % Concentration over time with first-order elimination

% Renal Clearance Profile
C_renal = C_p_initial * exp(-k_e_free * t); % Renal clearance with free drug assumption

% Hepatic Clearance Profile
C_hepatic = C_p_initial * exp(-Q_free / V_d * t); % Hepatic clearance with free drug assumption

% Michaelis-Menten Clearance Profile
% Adjust initial concentration for free drug
C_p_initial_free = C_p_initial * free_fraction;
[t_michaelis, C_michaelis] = ode45(@(t, Cp) -V_max * Cp / (K_m + Cp), t, C_p_initial_free);

% Plotting
figure;
hold on;

% Define HSL colors
colors = [hsl2rgb(180, 0.6, 0.6); % First Order Clearance (light blue)
          hsl2rgb(120, 0.6, 0.6); % Renal Clearance (light green)
          hsl2rgb(0, 0.6, 0.6);   % Hepatic Clearance (light red)
          hsl2rgb(240, 0.6, 0.6)]; % Michaelis-Menten Clearance (light purple)

% Plot each profile with the custom colors
%plot(t, C_first_order, 'Color', colors(1,:), 'LineWidth', 2, 'DisplayName', 'First Order Clearance');
plot(t, C_renal, 'Color', colors(2,:), 'LineWidth', 2, 'DisplayName', 'Renal Clearance');
plot(t, C_hepatic, 'Color', colors(3,:), 'LineWidth', 2, 'DisplayName', 'Hepatic Clearance');
%plot(t_michaelis, C_michaelis, 'Color', colors(4,:), 'LineWidth', 2, 'DisplayName', 'Michaelis-Menten Clearance');

hold off;

% Formatting the plot
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'LineWidth', 1.5);
xlabel('Time (hours)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Concentration (mg/L)', 'FontSize', 14, 'FontWeight', 'bold');
title('Effect of Plasma Protein Binding on Drug Clearance Mechanisms', 'FontSize', 16, 'FontWeight', 'bold');
legend('show', 'Location', 'northeast');
grid on;
xlim([0 24]);
ylim([0 max([C_first_order, C_renal, C_hepatic, C_michaelis]) * 1.1]); % Adjust y-axis to fit all profiles

% Function to convert HSL to RGB
function rgb = hsl2rgb(h, s, l)
    c = (1 - abs(2 * l - 1)) * s;
    x = c * (1 - abs(mod(h / 60, 2) - 1));
    m = l - c / 2;
    switch floor(h / 60)
        case 0
            rgb = [c, x, 0] + m;
        case 1
            rgb = [x, c, 0] + m;
        case 2
            rgb = [0, c, x] + m;
        case 3
            rgb = [0, x, c] + m;
        case 4
            rgb = [x, 0, c] + m;
        case 5
            rgb = [c, 0, x] + m;
    end
end

