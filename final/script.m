% Define the clearance rates for peripheral tissue and liver
k_pt = 0.05;  % 5% clearance from peripheral tissue each time step
k_liver = 0.10;  % 10% clearance from liver each time step

% Transition matrix with the new virtual clearance compartment
%     S, I, P, VLDL, LDL, HDL, L, Cl
A = [0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;  % Stomach transitions
     0.0, 0.5, 0.4, 0.05, 0.05, 0.0, 0.0, 0.0;  % Intestine transitions
     0.0, 0.0, 0.3, 0.1, 0.3, 0.2, 0.0, k_pt;  % Peripheral tissue transitions (clearing to clearance)
     0.0, 0.0, 0.0, 0.1, 0.8, 0.1, 0.0, 0.0;  % VLDL transitions
     0.0, 0.0, 0.0, 0.0, 0.2, 0.5, 0.3, 0.0;  % LDL transitions
     0.0, 0.0, 0.0, 0.0, 0.1, 0.7, 0.2, 0.0;  % HDL transitions
     0.05, 0.05, 0.05, 0.3, 0.3, 0.1, 0.1, k_liver;  % Liver transitions (clearing to clearance)
     0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.0ljj];  % Virtual clearance compartment


% Initial state vector (start with food in stomach)
x0 = zeros(length(A), 1);
x0(1) = 1e6; % Cholesterol/TG in stomach

% Simulate over time (10 time steps)
num_steps = 1000;
x = zeros(length(A), num_steps);  % Store state probabilities
x(:,1) = x0;  % Set initial state
dt = 0.5;
for t = 2:num_steps
    x(:, t) = A * x(:, t-1);  % Update state probabilities
end

% Plot state probabilities for each compartment
figure(1); clf;
plot(1:num_steps, x(8, :))
% plot(1:num_steps, x(1,:), 'r', 'LineWidth', 2); hold on;
% plot(1:num_steps, x(2,:), 'g', 'LineWidth', 2);
% plot(1:num_steps, x(3,:), 'b', 'LineWidth', 2);
% plot(1:num_steps, x(4,:), 'm', 'LineWidth', 2);
% plot(1:num_steps, x(5,:), 'c', 'LineWidth', 2);
% plot(1:num_steps, x(6,:), 'k', 'LineWidth', 2);
% plot(1:num_steps, x(7,:), 'y', 'LineWidth', 2);
legend('Stomach', 'Intestine', 'Peripheral Tissue', 'VLDL', 'LDL', 'HDL', 'Liver');
title('Cholesterol Distribution Over Time');
xlabel('Time Step');
ylabel('Probability');