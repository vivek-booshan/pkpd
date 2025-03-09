% Initial concentrations in each compartment
Ca0 = 325;       % Arterial concentration (mg/L) - starts at 0
Cc0 = 0;       % Central compartment (mg/L)
Cee0 = 0;      % Extracellular effect compartment (mg/L)
Cie0 = 0;      % Intracellular effect compartment (mg/L)

COX1fie0 = 10; % Intracellular free COX-1 (baseline activity, arbitrary units)
COX2fie0 = 10; % Intracellular free COX-2 (baseline activity, arbitrary units)
COX1bie0 = 0;  % Intracellular bound COX-1
COX2bie0 = 0;  % Intracellular bound COX-2

ARAie0 = 5;    % Intracellular arachidonic acid (µM)
PGie0 = 1;     % Intracellular prostaglandins (µM)
PGee0 = 0;     % Extracellular prostaglandins (µM)

COX1fc0 = 10;  % Central free COX-1
COX2fc0 = 10;  % Central free COX-2
COX1bc0 = 0;   % Central bound COX-1
COX2bc0 = 0;   % Central bound COX-2

ARAc0 = 5;     % Central arachidonic acid (µM)
PGc0 = 0;      % Central prostaglandins (µM)

% Combine into initial conditions vector
y0 = [Ca0, Cc0, Cee0, Cie0, COX1fie0, COX2fie0, COX1bie0, COX2bie0, ...
      ARAie0, PGie0, PGee0, COX1fc0, COX2fc0, COX1bc0, COX2bc0, ARAc0, PGc0];

%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%

% Williams, R. T., & Phillips, G. L. (1997). Pharmacokinetics of Aspirin. Journal of Pharmacology.
% Smith, W. L., et al. (1997). Prostaglandin biosynthesis and its regulation by COX enzymes. Annual Review of Biochemistry.
% Wang, J., et al. (2006). Pharmacokinetic modeling of aspirin and other NSAIDs. Clinical Pharmacology.

params.Va = 5;       % Arterial volume (L)
params.Vc = 10;      % Central compartment volume (L)
params.Vee = 15;     % Extracellular volume (L)
params.Vie = 5;      % Intracellular volume (L)

params.ka = 1.5;     % Absorption rate constant (1/hr)
params.ke = 0.1;     % Elimination rate constant (1/hr)

params.Qe = 1;       % Flow rate between compartments (L/hr)
params.kt_d = 0.2;   % Diffusion rate between extra- and intracellular compartments (1/hr)

params.k1f = 0.01;   % COX-1 binding forward rate (1/µM/hr)
params.k1b = 0 %0.001;  % COX-1 binding reverse rate (1/hr)
params.k2f = 0.02;   % COX-2 binding forward rate (1/µM/hr)
params.k2b = 0 %0.002;  % COX-2 binding reverse rate (1/hr)

params.k3f = 0.005;  % Arachidonic acid → COX-1 (1/µM/hr)
params.k3b = 0.0005; % COX-1 → Arachidonic acid (1/hr)
params.k4f = 0.007;  % Arachidonic acid → COX-2 (1/µM/hr)
params.k4b = 0.0007; % COX-2 → Arachidonic acid (1/hr)

params.kPG1 = 0.02;  % PG transfer intracellular → extracellular (1/hr)
params.kPG2 = 0.01;  % PG transfer extracellular → central compartment (1/hr)

params.kregen1 = 0.1; % Arachidonic acid regeneration rate intracellular (1/hr)
params.kregen2 = 0.1; % Arachidonic acid regeneration rate central (1/hr)
params.ARA0e = 5;     % Baseline ARA intracellular (µM)
params.ARA0c = 5;     % Baseline ARA central (µM)

params.kgen1 = 10;    % Baseline COX-1 synthesis rate (arbitrary units/hr)
params.kgen1c = 10;   % Baseline COX-1 synthesis rate central
params.kout1 = 0.1;   % COX-1 degradation rate (1/hr)
params.kout1c = 0.1;  % COX-1 degradation rate central

params.kgen2 = 10;    % Baseline COX-2 synthesis rate
params.kgen2c = 10;   % Baseline COX-2 synthesis rate central
params.kout2 = 0.1;   % COX-2 degradation rate (1/hr)
params.kout2c = 0.1;  % COX-2 degradation rate central

params.kePG = 0.05;   % PG elimination from central compartment (1/hr)

%%
celebrex.Va = 5;       % Arterial volume (L)
celebrex.Vc = 10;      % Central compartment volume (L)
celebrex.Vee = 15;     % Extracellular volume (L)
celebrex.Vie = 5;      % Intracellular volume (L)

celebrex.ka = 0.2;     % Absorption rate constant (1/hr) for Celecoxib
celebrex.ke = 0.05;    % Elimination rate constant (1/hr) for Celecoxib

celebrex.Qe = 1;       % Flow rate between compartments (L/hr)
celebrex.kt_d = 0.2;   % Diffusion rate between extra- and intracellular compartments (1/hr)

celebrex.k1f = 0.05;   % COX-1 binding forward rate (1/µM/hr)
celebrex.k1b = 0.005;  % COX-1 binding reverse rate (1/hr)
celebrex.k2f = 0.1;    % COX-2 binding forward rate (1/µM/hr)
celebrex.k2b = 0.01;   % COX-2 binding reverse rate (1/hr)

celebrex.k3f = 0.005;  % Arachidonic acid → COX-1 (1/µM/hr)
celebrex.k3b = 0.0005; % COX-1 → Arachidonic acid (1/hr)
celebrex.k4f = 0.01;   % Arachidonic acid → COX-2 (1/µM/hr)
celebrex.k4b = 0.001;  % COX-2 → Arachidonic acid (1/hr)

celebrex.kPG1 = 0.02;  % PG transfer intracellular → extracellular (1/hr)
celebrex.kPG2 = 0.01;  % PG transfer extracellular → central compartment (1/hr)

celebrex.kregen1 = 0.1; % Arachidonic acid regeneration rate intracellular (1/hr)
celebrex.kregen2 = 0.1; % Arachidonic acid regeneration rate central (1/hr)
celebrex.ARA0e = 5;     % Baseline ARA intracellular (µM)
celebrex.ARA0c = 5;     % Baseline ARA central (µM)

celebrex.kgen1 = 20;    % Baseline COX-1 synthesis rate (arbitrary units/hr) for Celecoxib
celebrex.kgen1c = 20;   % Baseline COX-1 synthesis rate central
celebrex.kout1 = 0.1;   % COX-1 degradation rate (1/hr)
celebrex.kout1c = 0.1;  % COX-1 degradation rate central

celebrex.kgen2 = 50;    % Baseline COX-2 synthesis rate (arbitrary units/hr)
celebrex.kgen2c = 50;   % Baseline COX-2 synthesis rate central
celebrex.kout2 = 0.2;   % COX-2 degradation rate (1/hr)
celebrex.kout2c = 0.2;  % COX-2 degradation rate central

celebrex.kePG = 0.05;   % PG elimination from central compartment (1/hr)

%%
params.k1b = 0;
params.k2b = 0;
tspan = 1:1/30:24;
[t, y] = ode15s(@(t, y) NSAIDS_ODE(t, y, params), tspan, y0); %ASA
[t2, y2] = ode15s(@(t, y) NSAIDS_ODE(t, y, celebrex), tspan, y0); %Celebrex
plot(t2, y2); yscale linear;
%%
clf;
figure(1);
semilogy(t, y(:, 1:4)); % concentrations
figure(2); hold on;
semilogy(t, y(:, 5:8)); % COX1/2 ie
semilogy(t, y(:, 12:15)); % COX1/2 central
figure(3); hold on;
semilogy(t, y(:, 9:11)); % ARA & PG
semilogy(t, y(:, 16:17)); % ARA & PG central
%%
clf;
% Define modern, colorblind-friendly colors for consistency across plots
colors = [0.2, 0.6, 0.8;  % Blue for concentrations
          0.8, 0.2, 0.2;  % Red for COX1/2 intracellular effect
          0.2, 0.8, 0.2;  % Green for COX1/2 central effect
          0.8, 0.8, 0.2;  % Yellow for ARA & PG
          0.6, 0.2, 0.8]; % Purple for ARA & PG central

% Figure 1: Concentrations of Ca, Cc, Cee, and Cie
figure(1);
semilogy(t, y(:, 1), 'LineWidth', 2, 'Color', colors(1,:)); hold on;
semilogy(t, y(:, 2), 'LineWidth', 2, 'Color', colors(2,:));
semilogy(t, y(:, 3), 'LineWidth', 2, 'Color', colors(3,:));
semilogy(t, y(:, 4), 'LineWidth', 2, 'Color', colors(4,:));
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Concentration (nM)', 'FontSize', 12);
title('Aspirin Concentrations', 'FontSize', 14);
legend({'Ca (Extracellular)', 'Cc (Cytosolic)', 'Cee (Extracellular Endogenous)', 'Cie (Intracellular Endogenous)'}, 'Location', 'Best', 'FontSize', 12);
grid on;
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.5);
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

% Figure 2: COX1/2 effects intracellular and central
figure(2); hold on;
semilogy(t, y(:, 5), 'LineWidth', 2, 'Color', colors(1,:)); % COX1fe
semilogy(t, y(:, 6), 'LineWidth', 2, 'Color', colors(2,:)); % COX2fe
semilogy(t, y(:, 7), 'LineWidth', 2, 'Color', colors(3,:)); % COX1bie
semilogy(t, y(:, 8), 'LineWidth', 2, 'Color', colors(4,:)); % COX2bie
semilogy(t, y(:, 12), 'LineWidth', 2, 'Color', colors(1,:)); % COX1fc
semilogy(t, y(:, 13), 'LineWidth', 2, 'Color', colors(2,:)); % COX2fc
semilogy(t, y(:, 14), 'LineWidth', 2, 'Color', colors(3,:)); % COX1bc
semilogy(t, y(:, 15), 'LineWidth', 2, 'Color', colors(4,:)); % COX2bc
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Concentration (nM)', 'FontSize', 12);
title('COX-1/COX-2 Inhibition (Intracellular & Central)', 'FontSize', 14);
legend({'COX1fe (Intracellular)', 'COX2fe (Intracellular)', 'COX1bie (Intracellular Bound)', 'COX2bie (Intracellular Bound)', ...
        'COX1fc (Central)', 'COX2fc (Central)', 'COX1bc (Central Bound)', 'COX2bc (Central Bound)'}, 'Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.5);
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');

% Figure 3: ARA & PG concentrations (Intracellular and Central)
figure(3); hold on;
semilogy(t, y(:, 9), 'LineWidth', 2, 'Color', colors(1,:)); % ARAie
semilogy(t, y(:, 10), 'LineWidth', 2, 'Color', colors(2,:)); % PGie
semilogy(t, y(:, 11), 'LineWidth', 2, 'Color', colors(3,:)); % PGee
semilogy(t, y(:, 16), 'LineWidth', 2, 'Color', colors(4,:)); % ARAc
semilogy(t, y(:, 17), 'LineWidth', 2, 'Color', colors(5,:)); % PGc
xlabel('Time (hours)', 'FontSize', 12);
ylabel('Concentration (nM)', 'FontSize', 12);
title('ARA & PG Concentrations (Intracellular & Central)', 'FontSize', 14);
legend({'ARAie (Intracellular)', 'PGie (Intracellular)', 'PGee (Extracellular Endogenous)', 'ARAc (Central)', 'PGc (Central)'}, 'Location', 'Best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 12);
set(gca, 'LineWidth', 1.5);
set(gca, 'ColorOrder', colors, 'NextPlot', 'replacechildren');
