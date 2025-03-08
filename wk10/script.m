clf; 

p.CL = 25.3; % liter/hr
p.Vc = 7.77; % liter
p.VP1 = 6.92; % liter
p.VP2 = 8.08; % liter
p.QP1 = 19.3; % liter/hr
p.QP2 = 2.33; % liter/hr
p.kgrowth = 1.46; % 1/hr
p.kdeath = 0.187; % 1/hr
p.Bmax = 5e8; % CFU/mL
p.Emax = 2.70; %1/hr
p.EC50 = 0.00531; %mg/liter
p.gamma = 1.06;


y0 = [0; 0; 0; 1e6; 0; 0];

[tcontrol, ycontrol] = ode45(@(t, y) penicillinODE(t, y, p), 0:1/60:24, y0);
[time, concentrations] = dothing(1000, 6, p);


%%
figure(1);
subplot(2, 1, 1);
plot(time, concentrations(:, 1), '-b', 'DisplayName', 'Cc (Blood Plasma)');
hold on;
plot(time, concentrations(:, 2), '-r', 'DisplayName', 'Cp1 (W.P.O.T)');
plot(time, concentrations(:, 3), '-g', 'DisplayName', 'Cp2 (P.P.T)');
xlabel('Time (hours)');
ylabel('Concentration (mg/L)');
legend();
title('Drug Concentrations Over Time');
grid on; yscale log;

subplot(2, 1, 2);
semilogy(time, concentrations(:, 4), '-', 'DisplayName', 'Bs (Drug Sensitive Bacteria)');
hold on;
semilogy(time, concentrations(:, 5), '-', 'DisplayName', 'Br (Drug Insensitive Bacteria)');
semilogy(tcontrol, ycontrol(:, 4), 'DisplayName', 'Control Bs');
semilogy(tcontrol, ycontrol(:, 5), 'DisplayName', 'Control Br');
xlabel('Time (hours)');
ylabel('Bacterial Population');
legend();
title('Bacterial Population Over Time');
grid on;
%%

[t6000, c6000] = dothing(6000, 1, p);
[t750, c750] = dothing(750, 8, p);
[t500, c500] = dothing(500, 12, p);
[t250, c250] = dothing(250, 24, p);

figure(2); clf; hold on;
plot(tcontrol, ycontrol(:, 4), 'DisplayName', 'Control Bs');
plot(t6000, c6000(:, 4), 'LineWidth', 1.5, 'DisplayName', '6000 mg (S)');
plot(time, concentrations(:, 4), 'LineWidth', 1.5, 'DisplayName', '1000 mg (S)');
plot(t750, c750(:, 4), 'LineWidth', 1.5, 'DisplayName', '750 mg (S)');
plot(t500, c500(:, 4), 'LineWidth', 1.5, 'DisplayName', '500 mg (S)');
plot(t250, c250(:, 4), 'LineWidth', 1.5, 'DisplayName', '250 mg (S)');

% Plot the drug-insensitive bacteria (R)
plot(tcontrol, ycontrol(:, 5), 'DisplayName', 'Control Br');
plot(t6000, c6000(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '6000 mg (R)');
plot(time, concentrations(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '1000 mg (R)');
plot(t750, c750(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '750 mg (R)');
plot(t500, c500(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '500 mg (R)');
plot(t250, c250(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '250 mg (R)');

% Add labels and legend
xlabel('Time (hours)');
ylabel('Bacteria Population (CFU/mL)');
title('Comparison of Bacterial Populations for Different Dosing Schedules');
legend('Location', 'best');
grid on; yscale log; ylim([1e0, max(ycontrol(:, 5))]);
hold off;

%%
figure(3); clf; hold on;

[t750, c750] = dothing(750, 6, p);
[t500, c500] = dothing(500, 6, p);
[t250, c250] = dothing(250, 6, p);

plot(tcontrol, ycontrol(:, 4), 'DisplayName', 'Control Bs');
plot(time, concentrations(:, 4), 'LineWidth', 1.5, 'DisplayName', '1000 mg (S)');
plot(t750, c750(:, 4), 'LineWidth', 1.5, 'DisplayName', '750 mg (S)');
plot(t500, c500(:, 4), 'LineWidth', 1.5, 'DisplayName', '500 mg (S)');
plot(t250, c250(:, 4), 'LineWidth', 1.5, 'DisplayName', '250 mg (S)');

% Plot the drug-insensitive bacteria (R)
plot(tcontrol, ycontrol(:, 5), 'DisplayName', 'Control Br');
plot(time, concentrations(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '1000 mg (R)');
plot(t750, c750(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '750 mg (R)');
plot(t500, c500(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '500 mg (R)');
plot(t250, c250(:, 5), '--', 'LineWidth', 1.5, 'DisplayName', '250 mg (R)');

xlabel('Time (hours)');
ylabel('Bacteria Population (CFU/mL)');
title('Comparison of Bacterial Populations for Different Doses');
legend('Location', 'best');
grid on; yscale log; ylim([1e0, max(ycontrol(:, 5))])
hold off;
%%
function [time, concentrations] = dothing(dose_mg, num_doses, p)

    total_time = 24;
    dose_interval = total_time / num_doses;
    
    dose_conc = dose_mg / p.Vc;
    
    y0 = [0; 0; 0; 1e6; 0; 0];
    
    tspan = 0:1/60:dose_interval;
    
    
    
    time = [];
    concentrations = [];
    %balance = [];
    
    for i = 1:num_doses
        y0(1) = y0(1) + dose_conc;
        [t, y] = ode45(@(t, y) penicillinODE(t, y, p), tspan, y0);
    
        time = [time(1:end-1); t + (i-1)*dose_interval];
        concentrations = [concentrations(1:end-1, :); y];
        %balance = [balance; dose_mg * i - p.Vc*y(:, 1) - p.VP1 * y(:,2) - p.VP2 * y(:, 3) - y(:, 6)];
        y0 = y(end, :);
    end
end