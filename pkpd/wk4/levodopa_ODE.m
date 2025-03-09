function dydt = levodopa_ODE(t, y, params)
    % Unpack parameters
    KTR_LCIG = params.KTR_LCIG;    % Transit absorption rate constant for LCIG (1/hr)
    KTR_LC_oral = params.KTR_LC_oral;  % Transit absorption rate constant for LC-oral (1/hr)
    CL_F = params.CL_F;            % Clearance (L/hr)
    Vc_F = params.Vc_F;            % Volume of central compartment (L)
    Q_F = params.Q_F;              % Inter-compartmental clearance (L/hr)
    Vp_F = params.Vp_F;            % Volume of peripheral compartment (L)
    F_rel_LCIG = params.F_rel_LCIG;% Bioavailability for LCIG relative to LC-oral
    
    % States
    A_gut = y(1);  % Drug amount in the gut (transit compartment)
    A_central = y(2);  % Drug amount in the central compartment
    A_peripheral = y(3);  % Drug amount in the peripheral compartment
    
    % Absorption process from the gut (transit absorption compartment)
    if strcmp(params.formulation, 'LCIG')
        dA_gut = -KTR_LCIG * A_gut;
        absorption_rate = KTR_LCIG * A_gut;
    else
        dA_gut = -KTR_LC_oral * A_gut;
        absorption_rate = KTR_LC_oral * A_gut;
    end
    
    % Central compartment (plasma)
    dA_central = absorption_rate * F_rel_LCIG - (CL_F / Vc_F) * A_central - (Q_F / Vc_F) * A_central + (Q_F / Vp_F) * A_peripheral;
    
    % Peripheral compartment (tissue distribution)
    dA_peripheral = (Q_F / Vc_F) * A_central - (Q_F / Vp_F) * A_peripheral;
    
    dA_elim = (CL_F/Vc_F)*A_central + (1-F_rel_LCIG)*absorption_rate;
    % Return derivatives
    dydt = [dA_gut; dA_central; dA_peripheral; dA_elim];
end
