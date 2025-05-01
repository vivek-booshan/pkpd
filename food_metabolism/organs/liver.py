import numpy as np
from .parameters import *
from .index import Index

def liver(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
    """
    The fat function computes the rate of change (dydt) for various metabolites 
    and processes in the adipose tissue (fat) compartment of a metabolic model.

    This function takes in the current state vector, `y`, representing concentrations or amounts 
    of metabolites, and updates it with the rate of change (`dydt`) based on metabolic interactions 
    governed by the parameters `p`.

    Attributes:
        t (float): The current time point (in the simulation's time units).
        y (np.ndarray): The state vector representing the current concentrations or amounts of metabolites.
        p (object): A parameters object containing rate constants and volume information.

    Returns:
        dydt (np.ndarray): The array of rate-of-change values for each state variable, which is used in 
                            numerical integration methods (e.g., `solve_ivp`) to simulate the system.
    """
    dydt = np.zeros(len(Index))
    __liver(t, y, p, dydt)
    return dydt

def __liver(t, y, p, dydt):
    """
    The fat function computes the rate of change (dydt) for various metabolites 
    and processes in the adipose tissue (fat) compartment of a metabolic model.

    This function takes in the current state vector, `y`, representing concentrations or amounts 
    of metabolites, and updates it with the rate of change (`dydt`) based on metabolic interactions 
    governed by the parameters `p`.

    Attributes:
        t (float): The current time point (in the simulation's time units).
        y (np.ndarray): The state vector representing the current concentrations or amounts of metabolites.
        p (object): A parameters object containing rate constants and volume information.
        dydt (np.ndarray): The array of rate-of-change values for each state variable, which is used in 
                            numerical integration methods (e.g., `solve_ivp`) to simulate the system.
    """
    __glucose(t, y, p, dydt)
    #__insulin(t, y, p, dydt)
    __fattyacids(t, y, p, dydt)
    __aminoacids(t, y, p, dydt)
    __g6p(t, y, p, dydt)
    __triglycerides(t, y, p, dydt)
    __pyruvate(t, y, p, dydt)
    __acetylcoa(t, y, p, dydt)
    __ROS(t, y, p, dydt)
    __fructose(t,y,p,dydt)
    return

def __glucose(t, y, p, dydt):
    dydt[Index.plasma_glucose] += (
        + (-p.Liver.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.Liver.k_G_to_plasma * y[Index.liver_glucose] * p.V.liver) / p.V.plasma
    )
    dydt[Index.liver_glucose] += (
        + (p.Liver.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.Liver.k_G_to_plasma * y[Index.liver_glucose] * p.V.liver) / p.V.liver
        - p.Liver.k_G_to_G6P * y[Index.liver_glucose] + p.Liver.k_G6P_to_G * y[Index.liver_G6P]
    )
    return

def __fructose(t,y,p,dydt):
    dydt[Index.plasma_fructose] +=(
        + (-p.Liver.k_F_from_plasma * y[Index.plasma_fructose] * p.V.plasma + p.Liver.k_F_to_plasma * y[Index.liver_fructose] * p.V.liver) / p.V.plasma
    )
    dydt[Index.liver_fructose]+=(
        + (p.Liver.k_F_from_plasma * y[Index.plasma_fructose] * p.V.plasma - p.Liver.k_F_to_plasma * y[Index.liver_fructose] * p.V.liver) / p.V.liver
        - p.Liver.k_F_to_P * y[Index.liver_fructose] 
    )

#def __insulin(t, y, p, dydt):
#    dydt[Index.plasma_insulin] += (
#        + (
#            - p.Liver.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma
#            + p.Liver.k_insulin_to_plasma * y[Index.liver_insulin] * p.V.liver
#        ) / p.V.plasma
#    )
#    
#    dydt[Index.liver_insulin] += (
#        + (p.Liver.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma - p.Liver.k_insulin_to_plasma * y[Index.liver_insulin] * p.V.liver) / p.V.liver
#        - p.Liver.kCL_insulin * y[Index.liver_insulin]
#    )
#    return

def __fattyacids(t, y, p, dydt):
    Km = 1
    
    dydt[Index.plasma_fattyacid] += (
        + (
            - p.Liver.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma
            + p.Liver.k_FA_to_plasma * y[Index.liver_fattyacid] * p.V.liver
        ) / p.V.plasma
    )

    dydt[Index.liver_fattyacid] += (
        + (p.Liver.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.Liver.k_FA_to_plasma * y[Index.liver_fattyacid] * p.V.liver) / p.V.liver
        - p.Liver.k_FA_to_ACoA * y[Index.liver_fattyacid]
        - 3 * (p.Liver.k_FA_to_TAG * p.V.liver * y[Index.liver_fattyacid] / (Km + y[Index.liver_fattyacid] * p.V.liver))**3
        + 3 * p.Liver.k_TAG_to_FA * y[Index.liver_TAG]
        + p.Liver.k_ACoA_to_FA * y[Index.liver_ACoA]
    )
    return

def __aminoacids(t, y, p, dydt):
    dydt[Index.plasma_aminoacid] += (
        + (
            - p.Liver.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            + p.Liver.k_AA_to_plasma * y[Index.liver_aminoacid] * p.V.liver
        ) / p.V.plasma
    )

    dydt[Index.liver_aminoacid] += (
        + (
            + p.Liver.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            - p.Liver.k_AA_to_plasma * y[Index.liver_aminoacid] * p.V.liver
        ) / p.V.liver
        - p.Liver.k_AA_to_ACoA * y[Index.liver_aminoacid]
    )
    return


def __g6p(t, y, p, dydt):
    dydt[Index.liver_G6P] += (
        + p.Liver.k_G_to_G6P * y[Index.liver_glucose]
        - p.Liver.k_G6P_to_G * y[Index.liver_G6P]
        - p.Liver.k_G6P_to_P * y[Index.liver_G6P]
        + p.Liver.k_P_to_G6P * y[Index.liver_pyruvate]**2
    )


def __triglycerides(t, y, p, dydt):
    Km = 1

    dydt[Index.liver_TAG] += (
        + (p.Liver.k_FA_to_TAG * p.V.liver * y[Index.liver_fattyacid] / (Km + y[Index.liver_fattyacid] * p.V.liver))**3
        - p.Liver.k_TAG_to_FA * y[Index.liver_TAG]
    )
    return


def __pyruvate(t, y, p, dydt):
    dydt[Index.liver_pyruvate] += (
        + 2 * p.Liver.k_G6P_to_P * y[Index.liver_G6P]
        - 2 * p.Liver.k_P_to_G6P * y[Index.liver_pyruvate]**2
        - p.Liver.k_P_to_ACoA * y[Index.liver_pyruvate]
        + p.Liver.k_F_to_P * y[Index.liver_fructose] 
    )
    return



def __acetylcoa(t, y, p, dydt):
    dydt[Index.liver_ACoA] += (
        + p.Liver.k_P_to_ACoA * y[Index.liver_pyruvate]
        + 8 * p.Liver.k_FA_to_ACoA * y[Index.liver_fattyacid]
        + p.Liver.k_AA_to_ACoA * y[Index.liver_aminoacid]
        - 8 * p.Liver.k_ACoA_to_FA * y[Index.liver_ACoA]
    )

def __ROS(t, y, p, dydt):
    # Fraction of reactions resulting in ROS (1% defa__ult)
    ROSpercent = 0.01
    Km = 1

    dydt[Index.liver_ROS] += ROSpercent * (
        p.Liver.k_FA_to_ACoA * y[Index.liver_fattyacid]
        + p.Liver.k_AA_to_ACoA * y[Index.liver_aminoacid]
        + 3 * (p.Liver.k_FA_to_TAG * p.V.liver * y[Index.liver_fattyacid] / (Km + y[Index.liver_fattyacid] * p.V.liver))**3
        + 3 * p.Liver.k_TAG_to_FA * y[Index.liver_TAG]
        + 8 * p.Liver.k_ACoA_to_FA * y[Index.liver_ACoA]
        + p.Liver.k_F_to_P * y[Index.liver_fructose] * 10
    )
    return
