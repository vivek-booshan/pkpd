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
p.krel = 0.6; % 1/hr

% Anatomical and Physiological Parameters for Tumor Tissue
p.Vcap = 0.0285; % (ml/g)
p.Vint = 0.323; % (ml/g)
p.Vecs = p.Vcap + p.Vint;
p.Vtu = 0.648; % (ml/g)
p.Q = 9.6; % (ml/hr/g)
p.k = 16.2; % [h-1/(ug/ml)]
p.ks = 0.0578; % 1/hr
p.fb = 0.2;
y0 = zeros(9, 1);
y0(5) = 15;
y0(9) = 1e10;
tspan = 1:1/60:48;
[time, soln] = ode45(@(t, y) DOX(t, y, p), tspan, y0);
semilogy(time, soln(:, [1, 3, 4, 5, 7]));
plot(time, soln(:, 9));
set(gca, 'YScale', 'log')
legend("Cb", "Cecs", "Ctu", "Cblipo", "Cintlipo");
title("Free and Liposomal DOX Concentrations")
xlabel("time (hrs)")
ylabel("Concentration (ml/g)")

%%
y0(5) = 11;
tspan = 1:1/60:200;
cell_count = [];
for i = 1:5
    y0(5) = y05(i);
    [time, soln] = ode45(@(t, y) DOX(t, y, p), tspan, y0);
    cell_count = [cell_count, soln(:, 9)];
end
%%
semilogy(time, cell_count)
title("Cell Count")
xlabel("time (hrs)")
ylabel("# of Cells")
legend("0 mg/kg", "6 mg/kg", "11 mg/kg", "15 mg/kg", "20 mg/kg");
ylim([1, 1e16])