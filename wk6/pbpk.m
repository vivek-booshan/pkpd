%% params

% Initialize volumes (Liters)
p.V.art = 1.5;       % Arterial blood volume
p.V.vein = 3.5;      % Venous blood volume
p.V.Gi = 1.5;        % Gastrointestinal volume
p.V.P = 0.1;         % Pancreas volume
p.V.Lvr = 1.5;       % Liver volume
p.V.Vsc = 2.0;       % Viscera volume
p.V.K = 0.3;         % Kidney volume
p.V.H = 0.3;         % Heart volume
p.V.B = 1.3;         % Brain volume
p.V.M = 30.0;        % Muscle volume
p.V.subq = 12.0;     % Subcutaneous tissue volume
p.V.L = 5.0;         % Lung volume

% Arterial to Organ Rate Constants (1/hr)
p.k_art.Gi = 1.5;    % Rate constant from arteries to GI
p.k_art.P = 1.2;     % Rate constant from arteries to pancreas
p.k_art.Lvr = 1.8;   % Rate constant from arteries to liver
p.k_art.Vsc = 1.6;   % Rate constant from arteries to viscera
p.k_art.K = 1.9;     % Rate constant from arteries to kidneys
p.k_art.H = 2.0;     % Rate constant from arteries to heart
p.k_art.B = 1.5;     % Rate constant from arteries to brain
p.k_art.M = 1.1;     % Rate constant from arteries to muscle
p.k_art.subq = 0.8;  % Rate constant from arteries to subcutaneous tissue
p.k_art.L = 1.5;     % Rate constant from arteries to lungs

% Venous Return Rate Constants (1/hr)
p.k_vein.Gi = 1.5;    % Rate constant from GI to veins
p.k_vein.P = 1.2;     % Rate constant from pancreas to veins
p.k_vein.Lvr = 1.8;   % Rate constant from liver to veins
p.k_vein.Vsc = 1.6;   % Rate constant from viscera to veins
p.k_vein.K = 1.9;     % Rate constant from kidneys to veins
p.k_vein.H = 2.0;     % Rate constant from heart to veins
p.k_vein.B = 1.5;     % Rate constant from brain to veins
p.k_vein.M = 1.1;     % Rate constant from muscle to veins
p.k_vein.subq = 0.8;  % Rate constant from subcutaneous tissue to veins
p.k_vein.L = 1.5;     % Rate constant from lungs to veins

% Diffusion btw Compartments (1/hr)
p.k_diff.M_subq = 0.1;
p.k_diff.subq_M = p.k_diff.M_subq / 100;

p.k_diff.Lvr_vsc = 0.15;
p.k_diff.vsc_Lvr = p.k_diff.Lvr_vsc;

p.k_diff.Gi_vsc = 0.15;
p.k_diff.vsc_Gi = p.k_diff.Gi_vsc / 10;

p.k_diff.vsc_P = 0.1;
p.k_diff.P_vsc = p.k_diff.vsc_P / 10;

p.k_diff.K_vsc = 0.01;
p.k_diff.vsc_K = p.k_diff.K_vsc / 10;

% Absorption Rate Constants (1/hr)
p.ka_Gi = 0.1;            % Absorption rate constant from GI
p.ka_M = 0;               % Absorption rate constant from muscle
p.ka_subq = 0;            % Absorption rate constant from subcutaneous tissue
p.ka_L = 0;               % Absorption rate constant from lungs

% Clearance Rate Constants (1/hr)
p.kcl_K = 0.12;       % Kidney clearance rate constant
p.kcl_Lvr = 0.15;     % Liver clearance rate constant
p.kcl_Gi = 0.05;      % GI clearance rate constant

% portal vein rates
p.klvr_Gi = 0.4;
p.klvr_P = 0.3;

% Display parameters
disp('PBPK model parameters:');
disp(p);
%%
y0 = zeros(1, 14);
dose = 10;
y0(14) = dose;
tspan = 1:1/60:24;

% Absorption Rate Constants (1/hr)
p.ka_Gi = 0.1;            % Absorption rate constant from GI
p.ka_M = 0;               % Absorption rate constant from muscle
p.ka_subq = 0;            % Absorption rate constant from subcutaneous tissue
p.ka_L = 0;               % Absorption rate constant from lungs

[t1, y1, b1] = solve_pbpk(@ode45, tspan, y0, p);
figure(1);
plot_pbpk(t1, y1, b1);

% Absorption Rate Constants (1/hr)
p.ka_Gi = 0;            % Absorption rate constant from GI
p.ka_M = 0.1;               % Absorption rate constant from muscle
p.ka_subq = 0;            % Absorption rate constant from subcutaneous tissue
p.ka_L = 0;               % Absorption rate constant from lungs

[t2, y2, b2] = solve_pbpk(@ode45, tspan, y0, p);
figure(2);
plot_pbpk(t2, y2, b2);

% Absorption Rate Constants (1/hr)
p.ka_Gi = 0;            % Absorption rate constant from GI
p.ka_M = 0;               % Absorption rate constant from muscle
p.ka_subq = 0.1;            % Absorption rate constant from subcutaneous tissue
p.ka_L = 0;               % Absorption rate constant from lungs

[t3, y3, b3] = solve_pbpk(@ode45, tspan, y0, p);
figure(3);
plot_pbpk(t3, y3, b3);

% Absorption Rate Constants (1/hr)
p.ka_Gi = 0;            % Absorption rate constant from GI
p.ka_M = 0;               % Absorption rate constant from muscle
p.ka_subq = 0;            % Absorption rate constant from subcutaneous tissue
p.ka_L = 0.1;               % Absorption rate constant from lungs

[t4, y4, b4] = solve_pbpk(@ode45, tspan, y0, p);
figure(4);
plot_pbpk(t4, y4, b4);
%%
figure(5);
tiledlayout(4,3, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% Labels for the absorption routes
absorptionRoutes = {'ka\_Gi', 'ka\_M', 'ka\_subq', 'ka\_L'};

% Labels for the compartments (assuming these are the names)
compartmentNames = {'Arteries', 'Gastrointestinal', 'Pancreas', 'Liver', 'Viscera', 'Kidney', 'Heart', 'Brain', 'Muscle', 'Subcutaneous', 'Lungs', 'Veins'};

% Loop through compartments and plot data
for i = 1:numel(y0) - 2
    nexttile; hold on;
    
    % Plot for each absorption route
    plot(t1, y1(:, i), 'LineWidth', 1.5);
    plot(t2, y2(:, i), 'LineWidth', 1.5);
    plot(t3, y3(:, i), 'LineWidth', 1.5);
    plot(t4, y4(:, i), 'LineWidth', 1.5);

    % Add labels for each subplot
    if (i == 10 || i == 11 || i == 12)
        xlabel('Time (hr)');
    end
    if mod(i, 3) == 1
        ylabel('Concentration (mmol/L)');
    end
    % Title for the compartment
    title(compartmentNames{i}, 'Interpreter', 'none');

    % Add the legend for the absorption routes
    legend(absorptionRoutes, 'Location', 'best');
    yscale log;
    xscale log;
    grid on;
    hold off;
end

sgtitle("Compartment")
%%
figure(1);
colors = [...
    0.8500, 0.3250, 0.0980;  % reddish-orange
    0.4660, 0.6740, 0.1880;  % green
    0.0000, 0.4470, 0.7410;  % blue
    0.4940, 0.1840, 0.5560;  % purple
    0.3010, 0.7450, 0.9330;  % light blue
    0.6350, 0.0780, 0.1840;  % dark red
    0.9290, 0.6940, 0.1250;  % yellowish-orange
    0.4660, 0.6740, 0.1880;  % olive green
    0.0000, 0.6196, 0.4500;  % teal green
    0.6902, 0.0196, 0.6902;  % magenta
    0.9608, 0.5098, 0.1882;  % dark orange
    0.3922, 0.8314, 0.0745
];

tiledlayout(4,1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
set(gcf, 'Position', [100, 100, 800, 600]); % Adjust the figure size

% Primary plot (3/4 of the figure) - Concentrations
nexttile([3, 1]);
hold on; 
for i = 1:size(colors, 1)
    plot(t, y(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
end
hold off;

xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
title('Concentrations in Organs Over Time');
legend({'Arteries', 'GI', 'Pancreas', 'Liver', 'Viscera', 'Kidney', 'Heart', 'Brain', 'Muscle', 'Subcutaneous', 'Lungs', 'Veins'}, ...
    'Location', 'eastoutside', 'FontSize', 9);
grid on; yscale log;
set(gca, 'FontSize', 12);

% Mass balance plot (1/4 of the figure) - Total Mass Over Time
nexttile([1, 1]);
plot(t, balance, 'k', 'LineWidth', 2); % Mass balance
xlabel('Time (hours)');
ylabel('Total Mass (mg)');
title('Mass Balance Over Time');
grid on;
set(gca, 'FontSize', 12);


