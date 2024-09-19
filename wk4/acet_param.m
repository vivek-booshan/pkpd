function p = acet_param(...
                            ka, kmg, kms, kmcm, ...
                            Vmg, Vms, Vmcm, V, Vg, Vs, ...
                            kep, keg, kes, ...
                            dose, Fep, tlag ...
                        )
    % create struct of parameters
    %
    % kmg, kms, kmcm (mmol/L) represent first order rate kinetics of
    % glucuronidation, sulphation, and cysteine/mercapture formation.
    %
    % vmg, vms, vmcm (mmol/h) represent Vmax for glucuronidation, sulphation and
    % cysteine/mercapture formation, respectively.
    %
    % v, vg, vs (L/kg) reprsent VoD of APAP, APAP-G, APAP-S
    %
    % kep, keg, kes (1/h) represent first order renal elimation of 
    % APAP, APAP-G, APAP-S.
    %
    % dose (mg/kg)
    % Fep fraction of paracetamol unaccounted for by metabolic pathways
    %
    % retrieved from 
    % https://onlinelibrary.wiley.com/doi/full/10.1111/j.1440-1681.2008.05029.x?saml_referrer=
    
    arguments
        ka   double = 2.04
        kmg  double = 7
        kms  double = 0.097
        kmcm double = 0.30
        Vmg  double = 64
        Vms  double = 740  
        Vmcm double = 0.25
        % V    double = 0.842 % L/kg 
        % Vg   double = 0.268 % L/kg
        % Vs   double = 17.5  % L/kg
        V    double = 56.4;
        Vg   double = 17.8;
        Vs   double = 17.5;
        kep  double = 0.07
        keg  double = 0.49
        kes  double = 2.44
        dose double = 60 %0.397
        Fep  double = 0.29
        tlag double = 0.54
    end

    p = struct('kmg', kmg, 'kms', kms, 'kmcm', kmcm, ...
        'Vmg', Vmg, 'Vms', Vms, 'Vmcm', Vmcm, ...
        'kep', kep, 'keg', keg, 'kes', kes, ...
        'V', V, 'Vg', Vg, 'Vs', Vs, ...
        'dose', dose, 'Fep', Fep, 'ka', ka, ...
        'tlag', tlag ...
        );
end