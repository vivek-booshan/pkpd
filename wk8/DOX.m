function dydt = DOX(t, y, p)
    Cb = y(1);
    X = y(2);
    Cecs = y(3);
    Ctu = y(4);
    Cblipo = y(5);
    Xres = y(6);
    Ccaplipo = y(7);
    Cintlipo = y(8);
    Cs = y(9); % Cell Concentrationw
    dydt = zeros(8, 1);

    % Free DOX
    % blood
    dydt(1) = (p.krel * p.Vclipo * Cblipo - p.k12 * p.Vc * Cb + p.k21 * X + p.Q * Cecs - p.Q * Cb)/p.Vc;
    % Tissue
    dydt(2) = p.k12 * p.Vc * Cb - p.k21 * X;
    % ECS
    dydt(3) = (p.Q*Cb + p.kte * p.Vtu * Ctu - p.Q * Cecs - p.ket * p.Vecs * Cecs + p.krel*(p.Vcap * Ccaplipo + p.Vint * Cintlipo))/p.Vecs;
    % Tumor
    dydt(4) = (p.ket * p.Vecs * Cecs - p.kte * p.Vtu * Ctu)/p.Vtu;
    

    % Liposomal DOX
    % Blood
    dydt(5) = (p.Q * Ccaplipo - p.Q * Cblipo - (p.kres + p.krel)*p.Vclipo * Cblipo)/p.Vclipo;
    % Res
    dydt(6) = p.kres * p.Vclipo * Cblipo;
    % CAP
    dydt(7) = (p.Q * Cblipo - p.Q * Ccaplipo - (p.krel + p.ktu) * p.Vcap * Ccaplipo) / p.Vcap;
    % Int
    dydt(8) = (p.ktu * p.Vcap * Ccaplipo - p.krel * p.Vint * Cintlipo)/p.Vint;
    
    dydt(9) = -p.k * p.fb * Cecs + p.ks * Cs;

end