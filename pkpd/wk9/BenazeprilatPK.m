function dydt = BenazeprilatPK(t, y, p)
    %
    % dydt(1) = Adepot
    % dydt(2) = Afree
    % dydt(3) = Abound
    %
    
    Tinf = 62/60;
    dydt = zeros(3, 1);
    if t <= Tinf
        dydt(1) = (p.D / Tinf) - p.ka * y(1);
    else
        dydt(1) = - p.ka * y(1);
    end
    dydt(2) = p.ka * y(1) - p.k10 * y(2) - p.k1 * y(2) * (p.BS - y(3)) + p.k2 * y(3);
    dydt(3) = p.k1 * y(2) * (p.BS - y(3)) - p.k2 * y(3);
end