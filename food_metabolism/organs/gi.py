import numpy as np
from organs.parameters import *
from organs.index import Index

def giInit():
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            subq=11.0,
            vsc=1.0,
            gut=1.25,
            liver=0.0,
            muscle=0.0,
            pancreas=0.0,
            brain=0.0
        ),
        Shared=SharedRates(
            k_P_to_ACoA=1.0,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=0.0,               # unused
            k_FA_to_ACoA=1 / 8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1 / 4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=1.0,              # acetylcoa_to_fattyacids
            k_G_to_G6P=1.0,                # glucose_to_g6p
            k_G6P_to_G=0.1,                # g6p_to_glucose
            k_P_to_G6P=0.1,                # pyruvate_to_g6p
            k_G6P_to_P=1.0                 # g6p_to_pyruvate
        ),
        Subq=None,
        Vsc=None,
        M=None,
        GI=GIParameters(
            kabs_glucose=0.1,
            kabs_fructose=0.1,
            kabs_fattyacid=0.05,
            k_diffusion_micelle_to_membrane=0.1,
            k_Vmax_trans=5.0,
            k_Vmax_reester=3.0,
            k_Vmax_export=2.0,
            Km_trans=10.0,
            Km_reester=5.0,
            Km_export=4.0,
            kCL_glucose=0.05,
            kCL_fructose=0.05,
            kCL_fattyacid=0.1
        )
        # Pancreas is not being used, so not included
    )
    return p

def GI(t, y, p):
    """
    The GI function computes the rate of change (dydt) for various metabolites in the gastrointestinal 
    (GI) compartment of a metabolic model, including glucose, fructose, and fatty acid metabolism.

    This function takes in the current state vector, `y`, representing concentrations or amounts of 
    metabolites, and updates it with the rate of change (`dydt`) based on metabolic interactions governed 
    by the parameters `p`.

    Attributes:
        t (float): The current time point (in the simulation's time units).
        y (np.ndarray): The state vector representing the current concentrations or amounts of metabolites in the GI tract.
        p (object): A parameters object containing rate constants and volume information.

    Returns:
        dydt (np.ndarray): The array of rate-of-change values for each state variable, used in numerical 
                            integration methods (e.g., `solve_ivp`) to simulate the system.
    """
    dydt = np.zeros(len(Index))
    glucose_two_compartment(t, y, p, dydt)
    fructose_two_compartment(t, y, p, dydt)
    fatty_acid_full_model(t, y, p, dydt)
    return dydt

def __GI(t, y, p, dydt):
    glucose_two_compartment(t, y, p, dydt)
    fructose_two_compartment(t, y, p, dydt)
    fatty_acid_full_model(t, y, p, dydt)

def glucose_two_compartment(t, y, p, dydt):
    V_gut = p.V.gut
    V_blood = p.V.plasma  # Using plasma as "blood" pool

    kabs = p.GI.kabs_glucose
    kclear = p.GI.kCL_glucose

    dydt[Index.gut_glucose] += -(kabs * y[Index.gut_glucose] * V_gut) / V_gut
    dydt[Index.plasma_glucose] += ((kabs * y[Index.gut_glucose] * V_gut) - (kclear * y[Index.plasma_glucose]* V_blood)) / V_blood

def fructose_two_compartment(t, y, p, dydt):
    V_gut = p.V.gut
    V_blood = p.V.plasma

    kabs = p.GI.kabs_fructose
    kclear = p.GI.kCL_fructose

    dydt[Index.gut_fructose] += -(kabs * y[Index.gut_fructose] * V_gut) / V_gut
    dydt[Index.plasma_fructose] += ((kabs * y[Index.gut_fructose] * V_gut) - (kclear * y[Index.plasma_fructose] * V_blood)) / V_blood

def fatty_acid_full_model(t, y, p, dydt):
    GI = p.GI
    V_blood = p.V.plasma

    # Pull out parameters
    J_diff = GI.k_diffusion_micelle_to_membrane * (y[Index.micellar_fattyacid] - y[Index.membrane_fattyacid])
    J_trans = (GI.k_Vmax_trans * y[Index.membrane_fattyacid]) / (GI.Km_trans + y[Index.membrane_fattyacid] + 1e-6)
    J_reester = (GI.k_Vmax_reester * y[Index.cytosol_fattyacid]) / (GI.Km_reester + y[Index.cytosol_fattyacid] + 1e-6)
    J_export = (GI.k_Vmax_export * y[Index.cytosol_TAG]) / (GI.Km_export + y[Index.cytosol_TAG] + 1e-6)
    J_clear = GI.kCL_fattyacid * y[Index.plasma_glucose]

    dydt[Index.micellar_fattyacid] += -J_diff
    dydt[Index.membrane_fattyacid] += J_diff - J_trans
    dydt[Index.cytosol_fattyacid] += J_trans - J_reester
    dydt[Index.cytosol_TAG] += J_reester - J_export
    dydt[Index.plasma_fattyacid] += (J_export - J_clear)
