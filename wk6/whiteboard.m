function [time, solution, balance] = solve_pbpk(odefun, tspan, y0, p)
    arguments
        odefun % integrator
        tspan (1, :) double
        y0 (1, :) double
        p struct
    end

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

    [time, solution] = odefun(@(t, y) pbpk_ODE(t, y, p), tspan, y0);
    balance = y0(14) - sum(V_vector' .* solution, 2);
end