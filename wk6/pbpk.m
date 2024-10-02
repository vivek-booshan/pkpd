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

% Absorption Rate Constants (1/hr)
p.ka_Gi = 0.1;        % Absorption rate constant from GI
p.ka_M = 0; %0.05;        % Absorption rate constant from muscle
p.ka_subq = 0; %0.05;     % Absorption rate constant from subcutaneous tissue
p.ka_L = 0; %0.1;         % Absorption rate constant from lungs

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
y0(14) = 2;

tspan = 1:1/60:24;
[t, y] = ode15s(@(t, y) pbpk_ODE(t, y, p), tspan, y0);
V_vector = [
    p.V.art;     % Arteries (index 1)
    p.V.Gi;     % Gastrointestinal (index 2)
    p.V.P;      % Pancreas (index 3)
    p.V.Lvr;    % Liver (index 4)
    p.V.Vsc;    % Viscera (index 5)
    p.V.K;      % Kidney (index 6)
    p.V.H;      % Heart (index 7)
    p.V.B;      % Brain (index 8)
    p.V.M;      % Muscle (index 9)
    p.V.subq;   % Subcutaneous (index 10)
    p.V.L;      % Lungs (index 11)
    p.V.vein;    % Veins (index 12)
    1;
    1;
];
balance = y0(14) - sum(V_vector'.*y, 2);
balance
plot(y(:, [1, 12]))
legend();

