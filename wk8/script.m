% Kinetic Parameters for Free DOX
p.Vc = 515; %(ml/250 g)
p.k12 = 5.13; % 1/hr
p.k21 = 0.927; % 1/hr
p.kel = 1.56; % 1/hr
p.ket = 2.27; % 1/hr
p.kte = 0.0485; % 1/hr

% Kinetic Parameters for Liposomal DOX
p.Vclipo = 16.7; % (ml/250g)
p.ktu = 7.17; % 1/hr
p.kres = 0.06; % 1/hr
p.krel = 0.06; % 1/hr

% Anatomical and Physiological Parameters for Tumor Tissue
p.Vcap = 0.0285; % (ml/g)
p.Vint = 0.323; % (ml/g)
p.Vecs = p.Vcap + p.Vint;
p.Vtu = 0.648; % (ml/g)
p.Q = 9.72; % (ml/hr/g)
p.k = 16.2; % [h-1/(ug/ml)]
p.ks = 0.0193; % 1/hr
p.fb = 0.2;
y0 = zeros(9, 1);
y0(1) = 6;
y0(5) = y0(1);
y0(9) = 100;
tspan = 1:1/60:60;
[time, soln] = ode45(@(t, y) DOX(t, y, p), tspan, y0);
%semilogy(time, soln(:, [1, 3, 4, 5, 7]));
plot(time, soln(:, 9));
legend("Cb", "Cecs", "Ctu", "Cblipo", "Cintlipo");