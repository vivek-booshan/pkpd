# TODO : need to include vldl incorporation?
import numpy as np
from .parameters import *
from .index import Index

def GI(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
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
    __GI(t, y, p, dydt)
    return dydt

def __GI(t, y, p, dydt):
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
        dydt (np.ndarray): The array of rate-of-change values for each state variable, used in numerical 
                            integration methods (e.g., `solve_ivp`) to simulate the system.
    """
    __glucose_two_compartment(t, y, p, dydt)
    __fructose_two_compartment(t, y, p, dydt)
    __fatty_acid_full_model(t, y, p, dydt)

def __glucose_two_compartment(t, y, p, dydt):
    V_gut = p.V.gut
    V_blood = p.V.plasma  # Using plasma as "blood" pool

    kabs = p.GI.kabs_glucose
    kclear = p.GI.kCL_glucose

    dydt[Index.gut_glucose] += -(kabs * y[Index.gut_glucose] * V_gut) / V_gut
    dydt[Index.plasma_glucose] += ((kabs * y[Index.gut_glucose] * V_gut) - (kclear * y[Index.plasma_glucose]* V_blood)) / V_blood

def __fructose_two_compartment(t, y, p, dydt):
    V_gut = p.V.gut
    V_blood = p.V.plasma

    kabs = p.GI.kabs_fructose
    kclear = p.GI.kCL_fructose

    dydt[Index.gut_fructose] += -(kabs * y[Index.gut_fructose] * V_gut) / V_gut
    dydt[Index.plasma_fructose] += ((kabs * y[Index.gut_fructose] * V_gut) - (kclear * y[Index.plasma_fructose] * V_blood)) / V_blood

def __fatty_acid_full_model(t, y, p, dydt):
    GI = p.GI
    V_blood = p.V.plasma

    # Pull out parameters
    J_diff = GI.k_diffusion_micelle_to_membrane * (y[Index.micellar_fattyacid] - y[Index.membrane_fattyacid])
    J_trans = (GI.k_Vmax_trans * y[Index.membrane_fattyacid]) / (GI.Km_trans + y[Index.membrane_fattyacid] + 1e-6)
    J_reester = (GI.k_Vmax_reester * y[Index.cytosol_fattyacid]) / (GI.Km_reester + y[Index.cytosol_fattyacid] + 1e-6)
    J_export = (GI.k_Vmax_export * y[Index.cytosol_TAG]) / (GI.Km_export + y[Index.cytosol_TAG] + 1e-6)
    J_clear = 0

    dydt[Index.micellar_fattyacid] += -J_diff
    dydt[Index.membrane_fattyacid] += J_diff - J_trans
    dydt[Index.cytosol_fattyacid] += J_trans - J_reester
    dydt[Index.cytosol_TAG] += J_reester - J_export
    dydt[Index.plasma_fattyacid] += J_trans
