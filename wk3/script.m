%% IV

%% Oral
q = 1; % mmol/hr
V = 1; % liters
kcl = 1; % 1/hr
ka = 1; % mmol/hr
oral = Model(q, V, kcl, ka);
y0 = [oral.params.q, 0, 0];

[t, y] = ode45(@Model.oralODE, [0, 10], y0, [], oral.params);

%[t, y] = Model.euler(@(t, y) Model.oralODE(t, y, oral.params), [0, 10], y0);

b = 1 - sum(y, 2);
plot(t, y)
% 
% analytic = @(t) ((q * ka)/(kcl + ka))*(exp(kcl*t) - exp(ka*t));
% plot(t, analytic(t))