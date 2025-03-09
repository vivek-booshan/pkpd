function dydt = penicillinODE(t, y, p)

    n=6;
    dydt = zeros(n, 1);
   
    %Cc
    dydt(1) = (- (p.CL / p.Vc) * y(1) +...
               - (p.QP1 / p.Vc) * y(1) + (p.QP1 / p.VP1) * y(2) + ...
               - (p.QP2 / p.Vc) * y(1) + (p.QP2 / p.VP2) * y(3));
    %Cp1
    dydt(2) = ((p.QP1 / p.Vc) * y(1) - (p.QP1 / p.VP1) * y(2));
    %Cp2
    dydt(3) = ((p.QP2 / p.Vc) * y(1) - (p.QP2 / p.VP2) * y(3));
    

    p.kSR = (y(4) + y(5))*(p.kgrowth - p.kdeath)/p.Bmax;
    DRUG = p.Emax * y(1)^p.gamma / (y(1)^p.gamma + p.EC50^p.gamma);

    % drug sensitive stage
    dydt(4) = p.kgrowth * y(4) - (p.kdeath + DRUG) * y(4) - p.kSR * y(4);
    
    % drug insensitive stage
    dydt(5) = p.kSR * y(4) - p.kdeath * y(5);
    
    % clearance
    dydt(6) = p.CL * y(1);
end