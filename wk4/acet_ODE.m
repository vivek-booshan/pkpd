% function dydt = acet_ODE(t, y, p)
%     % acetaminophen ode retrieved from 
%     % https://onlinelibrary.wiley.com/doi/full/10.1111/j.1440-1681.2008.05029.x?saml_referrer=
%     % where G, S, and CM stand for glucuronidation, sulphation, and
%     % cysteine/mercapture formation, respectively
%     %
%     % Inputs
%     %   t (float) : dummy input to log timestep
%     %   y (7x1 float) : (
%     %       plasma concentration
%     %                   1. Cp  ; APAP
%     %                   2. CpG ; APAP-G 
%     %                   3. CpS ; APAP-S
%     %       cumulative urine recovery
%     %                   4. UP  ; APAP
%     %                   5. UG  ; APAP-G
%     %                   6. US  ; APAP-S
%     %                   7. UCM ; APAP-CM
%     %   )
%     %   p (struct) : parameter struct (see acet_param.m for more)
% 
% 
% 
%     if p.tlag >= t 
%         input = p.ka * (p.dose/p.V) * exp(-p.ka*(t - p.tlag));
%     else
%         input = 0;
%     end
% 
%     %input = p.dose/p.V;
%     dydt = zeros(7, 1);
%     dydt(1) = input - y(1) * p.Vmg / (p.V * (p.kmg + y(1))) ...
%                     - y(1) * p.Vms / (p.V * (p.kms + y(1))) ...
%                     - y(1) * p.Vmcm / (p.V * (p.kmcm + y(1))) - y(1) * p.kep;
% 
%     dydt(2) =       + y(1) * p.Vmg / (p.Vg * (p.kmg + y(1))) - y(2) * p.keg;
%     dydt(3) =       + y(1) * p.Vms / (p.Vs * (p.kms + y(1))) - y(3) * p.kes;
% 
%     dydt(4) =       + y(1) * p.kep * p.V * p.Fep;
%     dydt(5) =       + y(2) * p.keg * p.Vg;
%     dydt(6) =       + y(3) * p.kes * p.Vs;
%     dydt(7) =       + y(1) * p.Vmcm / (p.kmcm + y(1));
% end

function dydt = acet_ODE(t, y, p)
    % Unpack parameters from structure
    KEP = p.KEP;
    KMG = p.KMG;
    VMG = p.VMG;
    KEG = p.KEG;
    KMS = p.KMS;
    VMS = p.VMS;
    KES = p.KES;
    VMCM = p.VMCM;
    FEP = p.FEP;
    KMCM = p.KMCM;
    DOSE1 = p.DOSE1;
    DOSE2 = p.DOSE2;
    TLAG1 = p.TLAG1;
    TLAG2 = p.TLAG2;
    V = p.V;
    KA = p.KA;

    % Compute inputs with time lag
    if t >= TLAG1
        INPUT1 = KA * (DOSE1 / V) * exp(-KA * (t - TLAG1));
    else
        INPUT1 = 0;
    end

    if t >= TLAG2
        INPUT2 = KA * (DOSE2 / V) * exp(-KA * (t - TLAG2));
    else
        INPUT2 = 0;
    end

    % Initialize the differential equations
    dydt = zeros(14, 1);

    % Differential equations
    dydt(1) = INPUT1 - y(1) * KEP - y(1) * VMG / ((KMG + y(1)) * V) - y(1) * VMS / ((KMS + y(1)) * V) - y(1) * VMCM / ((KMCM + y(1)) * V);
    dydt(2) = y(1) * VMG / ((KMG + y(1)) * (V * 0.26)) - y(2) * KEG;
    dydt(3) = y(1) * VMS / ((KMS + y(1)) * (V * 0.26)) - y(3) * KES;
    dydt(4) = INPUT2 - y(4) * KEP - y(4) * VMG / ((KMG + y(4)) * V) - y(4) * VMS / ((KMS + y(4)) * V) - y(4) * VMCM / ((KMCM + y(4)) * V);
    dydt(5) = y(4) * VMG / ((KMG + y(4)) * (V * 0.26)) - y(5) * KEG;
    dydt(6) = y(4) * VMS / ((KMS + y(4)) * (V * 0.26)) - y(6) * KES;
    dydt(7) = y(1) * KEP * FEP * V;
    dydt(8) = y(2) * KEG * (V * 0.26);
    dydt(9) = y(3) * KES * (V * 0.26);
    dydt(10) = y(1) * VMCM / (KMCM + y(1));
    dydt(11) = y(4) * KEP * FEP * V;
    dydt(12) = y(5) * KEG * (V * 0.26);
    dydt(13) = y(6) * KES * (V * 0.26);
    dydt(14) = y(4) * VMCM / (KMCM + y(4));
end




