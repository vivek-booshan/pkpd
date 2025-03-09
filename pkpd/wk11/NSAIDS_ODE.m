function dydt = NSAIDS_ODE(t, y, params)
    dydt = zeros(17, 1);

    Ca = y(1);        dCadt = dydt(1);
    Cc = y(2);        dCcdt = dydt(2);
    Cee = y(3);       dCeedt = dydt(3);
    Cie = y(4);       dCiedt = dydt(4);
    
    COX1fie = y(5);   dCOX1fiedt = dydt(5);
    COX2fie = y(6);   dCOX2fiedt = dydt(6);
    COX1bie = y(7);   dCOX1biedt = dydt(7);
    COX2bie = y(8);   dCOX2biedt = dydt(8);

    ARAie = y(9);     dARAiedt = dydt(9);
    PGie = y(10);     dPGiedt = dydt(10);
    PGee = y(11);     dPGeedt = dydt(11);

    COX1fc = y(12);   dCOX1fcdt = dydt(11);
    COX2fc = y(13);   dCOX2fcdt = dydt(12);
    COX1bc = y(14);   dCOX1bcdt = dydt(13);
    COX2bc = y(15);   dCOX2bcdt = dydt(14);
    
    ARAc = y(16);     dARAcdt = dydt(15);
    PGc = y(17);      dPGcdt = dydt(16);

    % Unpack parameters
    Va = params.Va; % Arterial volume
    Vc = params.Vc; % Central compartment volume
    Vee = params.Vee; % Extracellular interstitial volume
    Vie = params.Vie; % Intracellular interstitial volume
    ka = params.ka; % Absorption rate constant
    Qe = params.Qe; % Flow rate between central and extracellular interstitial compartments
    kt_d = params.kt_d; % Diffusion rate between extracellular and intracellular compartments
    k1f = params.k1f; % Forward binding rate for COX-1
    k1b = params.k1b; % Reverse binding rate for COX-1
    k2f = params.k2f; % Forward binding rate for COX-2
    k2b = params.k2b; % Reverse binding rate for COX-2
    k3f = params.k3f; % ARA --> COX1
    k3b = params.k3b; % COX1 --> ARA
    k4f = params.k4f; % ARA --> COX12
    k4b = params.k4b; % COX2 --> ARA
    ke = params.ke; % Elimination rate from the central compartment
    kgen1 = params.kgen1;
    kgen1c = params.kgen1c;
    kout1 = params.kout1;
    kout1c = params.kout1c;
    kgen2 = params.kgen2;
    kgen2c = params.kgen2c;
    kout2 = params.kout2;
    kout2c = params.kout2c;
    ARA0e = params.ARA0e;
    kPG1 = params.kPG1;
    kPG2 = params.kPG2;
    kregen1 = params.kregen1;
    ARA0c = params.ARA0c;
    kePG = params.kePG;

    dCadt = -ka * Ca / Va;
    dCcdt = ka * Va / Vc * Ca ...
        - Qe / Vc * (Cc - Cee) ...
        - ke * Cc ...
        - k1f * Cc * COX1fc + k1b * COX1bc ...
        - k2f * Cc * COX2fc + k2b * COX2bc;
    dCeedt = Qe / Vc * (Cc - Cee) ...
        - kt_d / Vee * (Cee - Cie);
    dCiedt = kt_d / Vie * (Cee - Cie) ...
        - k1f * Cie * COX1fie + k1b * COX1bie ...
        - k2f * Cie * COX2fie + k2b * COX2bie;

    dCOX1fiedt = kgen1 - kout1*(COX1fie / (COX1fie + COX1bie)) - k1f * Cie * COX1fie + k1b * COX1bie;
    dCOX2fiedt = kgen2 - kout2 * (COX2fie / (COX2fie + COX2bie)) - k2f * Cie * COX2fie + k2b * COX2bie;
    dCOX1biedt = -kout1 * (COX1bie / (COX1fie + COX1bie)) + k1f * Cie * COX1fie - k1b * COX1bie;
    dCOX2biedt = -kout2 * (COX2bie / (COX2fie + COX2bie)) + k2f * Cie * COX2fie - k2b * COX2bie;

    dARAiedt = kregen1 * (ARA0e - ARAie) - k3f * ARAie * COX1fie - k4f * ARAie * COX2fie + (k3b + k4b) * PGie;
    dPGiedt = k3f * ARAie * COX1fie + k4f * ARAie * COX2fie - (k3b + k4b)*PGie - kPG1 * (PGie - PGee);
    dPGeedt = kPG1 * (PGie - PGee) - kPG2 * (PGee - PGc);

    dCOX1fcdt = kgen1c - kout1c * (COX1fc / (COX1fc + COX1bc)) - k1f * Cc * COX1fc + k1b * COX1bc;
    dCOX2fcdt = kgen2c - kout2c * (COX2fc / (COX2fc + COX2bc)) - k2f * Cc * COX2fc + k2b * COX2bc;

    dCOX1bcdt = -kout1 * (COX1bc / (COX1fc + COX1bc)) + k1f * Cc * COX1fc - k1b * COX1bc;
    dCOX2bcdt = -kout2 * (COX2bc / (COX2fc + COX2bc)) + k2f * Cc * COX2fc - k1b * COX2bc;

    dARAcdt = kregen1 * (ARA0c - ARAc) - k3f * ARAc * COX1fc - k4f * ARAc * COX2fc + (k3b + k4b) * PGc;
    dPGcdt = k3f * ARAc * COX1fc + k4f * ARAc * COX2fc - (k3b + k4b) * PGc + kPG2  * (PGee - PGc) - kePG * PGc;

    dydt = [dCadt, dCcdt, dCeedt, dCiedt, dCOX1fiedt, dCOX2fiedt, dCOX1biedt, dCOX2biedt, ...
        dARAiedt, dPGiedt, dPGeedt, dCOX1fcdt, dCOX2fcdt, dCOX1bcdt, dCOX2bcdt, dARAcdt, dPGcdt]';
end