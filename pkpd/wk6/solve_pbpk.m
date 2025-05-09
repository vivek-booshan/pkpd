function [time, solution, balance] = solve_pbpk(odefun, tspan, y0, p)
    arguments
        odefun % integrator
        tspan (1, :) double
        y0 (1, :) double
        p struct
        %plot logical
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
    
    figure;
     colors = [...
        0.8500, 0.3250, 0.0980;  % reddish-orange
        0.4660, 0.6740, 0.1880;  % green
        0.0000, 0.4470, 0.7410;  % blue
        0.4940, 0.1840, 0.5560;  % purple
        0.9290, 0.6940, 0.1250;  % yellow
        0.6350, 0.0780, 0.1840;  % dark red
        0.3010, 0.7450, 0.9330;  % cyan
        0.6902, 0.0196, 0.6902;  % magenta
        0.9608, 0.5098, 0.1882;  % dark orange
        0.3922, 0.8314, 0.0745;  % lime green
        0.7039, 0.3216, 0.1569;  % brown
        0.9569, 0.7529, 0.3922   % peach
    ];
    
    tiledlayout(4,1, 'TileSpacing', 'Compact', 'Padding', 'Compact');
    set(gcf, 'Position', [100, 100, 800, 600]); % Adjust the figure size
    
    % Primary plot (3/4 of the figure) - Concentrations
    nexttile([3, 1]);
    hold on; 
    for i = 1:size(colors, 1)
        plot(time, solution(:, i), 'Color', colors(i, :), 'LineWidth', 1.5);
    end
    hold off;
    
    xlabel('Time (hours)');
    ylabel('Concentration (mg/L)');
    title('Concentrations in Organs Over Time');
    legend({'Arteries', 'GI', 'Pancreas', 'Liver', 'Viscera', 'Kidney', 'Heart', 'Brain', 'Muscle', 'Subcutaneous', 'Lungs', 'Veins'}, ...
        'Location', 'eastoutside', 'FontSize', 9);
    grid on; 
    yscale log;
    set(gca, 'FontSize', 12);
    
    % Mass balance plot (1/4 of the figure) - Total Mass Over Time
    nexttile([1, 1]);
    plot(time, balance, 'k', 'LineWidth', 2); % Mass balance
    xlabel('Time (hours)');
    ylabel('Total Mass (mg)');
    title('Mass Balance Over Time');
    grid on;
    set(gca, 'FontSize', 12);
end