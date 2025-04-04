def ligand_bound(L, R, kon, koff):
    """
    % Assumes Pseudo steady state assumption
    arguments (Input)
        L    double % Free Ligand Concentration
        R    double % Free Receptor Concentration
        kon  double % L + R --> LR
        koff double % LR --> L + R
    end
    arguments (Output)
        LR double % Bound Ligand Concentration
    end
    """
    LR = kon * L * R / koff
    return LR

def fractional_occupancy_Kd(L, Kd):
    f = L / (L + Kd)
    return f

def fractional_occupancy_LR(R, LR):
    f = 1 / (1 + (R / LR))
    return f

def fractional_inhibitory_effect(f_occupancy, k):
    denom = 1 + exp(k * (f - 0.5))
    f = 1/denom
    return f

def inhibitory_secretion(L, R, LR, kon, koff, kmax_secretion):
    Kd = L * R / LR
    fo = L / (L + Kd)
    f = 1 / (1 + exp(k * (f - 0.5)))
    keff = kmax_secretion * f
    return keff

def ligandODE(L, R, LR, kon, koff):
    dL = -kon * L * R + koff * LR
    dR = -kon * L * R + koff * LR
    dLR = +kon * L * R - koff * lR
    return (dL, dR, dLR)


