function out = BenazeprilatPD(delta, y, f)
    
    Imax_AII = 1;
    IC50_AII = 17.9;
    gamma_AII = 1;

    Emax_RA = 2.8;
    EC50_RA = 1;
    gamma_RA = 1.7;

    Imax_ALD = 0.9;
    IC50_ALD = 3.6;
    gamma_ALD = 2.7;

    % Angiotensin II
    f_AII = f(:, 1);
    f_RA = f(:, 2);
    f_ALD = f(:, 3);

    AII = y(:, 2) + y(:, 3);

    E1 = 1 - ((Imax_AII * AII.^gamma_AII)./(IC50_AII.^gamma_AII + AII.^gamma_AII));
    size(E1)
    out(:, 1) = f_AII .* E1;

    % Renin Activity
    E2 = 1 + ((Emax_RA * delta.^gamma_RA)./(EC50_RA^gamma_RA + delta.^gamma_RA));
    out(:, 2) = f_RA .* E2;

    % Aldosterone Data
    E3 = 1 - ((Imax_ALD * delta.^gamma_ALD)./(IC50_ALD^gamma_ALD + delta.^gamma_ALD));
    out(:, 3) = f_ALD .* E3;
end