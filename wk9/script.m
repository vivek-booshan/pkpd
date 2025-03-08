placebo_AII.M = 8.2;
placebo_AII.A = 2.4;
placebo_AII.phi = (23*60 + 20)/(24*60);

placebo_RA.M = 372;
placebo_RA.A = 108;
placebo_RA.phi = (23*60 + 20)/(24*60);

placebo_ALD.M = 245;
placebo_ALD.A = 71;
placebo_ALD.phi = (23*60 + 20)/(24*60);

tspan = 1:1/60:24;

p.ka = 2.2;
p.k10 = 0.56;
p.k1 = 0.011;
p.k2 = 0.010;
p.BS = 0.025;

y0 = zeros(3, 1);
p.D = 5; % 5 mg

[t, y] = ode45(@(t, y) BenazeprilatPK(t, y, p), tspan, y0);
delta = (t <= (62/60)); 
f(:, 1) = predictedPlaceboValue(t, placebo_AII);
f(:, 2) = predictedPlaceboValue(t, placebo_RA);
f(:, 3) = predictedPlaceboValue(t, placebo_ALD);

out = BenazeprilatPD(delta, y, f);