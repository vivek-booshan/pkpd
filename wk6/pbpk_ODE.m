function dydt = pbpk_ODE(t, y, p)

    dydt = zeros(14, 1);
    
    % Arteries
    name = fieldnames(p.k_art);
    sum_kart = 0;
    for i = 1:numel(name)
        sum_kart = sum_kart + p.V.(name{i});
    end
    dydt(1) = -p.V.art*y(1)*sum_kart + 2*p.k_art.L * p.V.L * y(11); % acct for subtraction of lung in sum so add double


    % Gastrointestinal (Gi)
    dydt(2) = p.ka_Gi * y(14);
    dydt(2) = dydt(2) + p.k_art.Gi * p.V.art * y(1) - p.klvr_Gi * p.V.Gi * y(2) - p.kcl_Gi * p.V.Gi * y(2);
    dydt(2) = dydt(2)/p.V.Gi;

    % Pancreas
    
    dydt(3) = p.k_art.P * p.V.art * y(1) - p.klvr_P * p.V.P * y(3);
    dydt(3) = dydt(3)/p.V.P;

    % Liver
    
    dydt(4) = p.k_art.Lvr * p.V.art * y(1) + p.klvr_Gi * p.V.Gi * y(2) + p.klvr_P * p.V.P * y(3) - p.k_vein.Lvr * p.V.Lvr * y(4) - p.kcl_Lvr * p.V.Lvr * y(4);
    dydt(4) = dydt(4)/p.V.Lvr;

    % Viscera
    
    dydt(5) = p.k_art.Vsc * p.V.art * y(1) - p.k_vein.Vsc * p.V.Vsc * y(5);
    dydt(5) = dydt(5)/p.V.Vsc;

    % Kidney
    
    dydt(6) = p.k_art.K * p.V.art * y(1) - p.k_vein.K * p.V.K * y(6) - p.kcl_K * p.V.K * y(6); 
    dydt(6) = dydt(6)/p.V.K;

    % Heart
    
    dydt(7) = p.k_art.H * p.V.art * y(1) - p.k_vein.H * p.V.H * y(7);
    dydt(7) = dydt(7)/p.V.H;

    % Brain
    
    dydt(8) = p.k_art.B * p.V.art * y(1) - p.k_vein.B * p.V.B * y(8);
    dydt(8) = dydt(8)/p.V.B;
    
    % Muscle

    dydt(9) = p.ka_M * y(14);
    dydt(9) = dydt(9) + p.k_art.M * p.V.art * y(1) - p.k_vein.M * p.V.M * y(9);
    dydt(9) = dydt(9)/p.V.M;
    
    % Subq
    dydt(10) = p.ka_subq * y(14);
    dydt(10) = dydt(10) + p.k_art.subq * p.V.art * y(1) - p.k_vein.subq * p.V.subq * y(10);
    dydt(10) = dydt(10)/p.V.subq;
    
    % Lungs

    dydt(11) = p.ka_L * y(14);
    dydt(11) = dydt(11) + p.k_vein.L * p.V.vein * y(1) - p.k_art.L * p.V.L * y(11);
    dydt(11) = dydt(11)/p.V.L;

    % Veins
    
    dydt(12) = p.k_vein.Gi * p.V.Gi * y(2) + p.k_vein.P * p.V.P * y(3) + p.k_vein.Lvr * p.V.Lvr * y(4) + ...
         p.k_vein.Vsc * p.V.Vsc * y(5) + p.k_vein.K * p.V.K * y(6) +  p.k_vein.H * p.V.H * y(7) + ...
          p.k_vein.B * p.V.B * y(8) + p.k_vein.M * p.V.M * y(9) + p.k_vein.subq * p.V.subq * y(10) + ...
          -p.k_vein.L * p.V.vein * y(11);

    % Clear (amt)

    dydt(13) = p.kcl_K * p.V.K * y(6) + p.kcl_Lvr * p.V.Lvr * y(4) + p.kcl_Gi * p.V.Gi * y(2);

    % Dose (amt)
    
    dydt(14) = -(p.ka_M + p.ka_Gi + p.ka_subq + p.ka_L) * y(14);

end