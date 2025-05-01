import numpy as np
from .parameters import *
from .index import Index

def fat(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
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
    __fat(t, y, p, dydt)
    return dydt

def __fat(t, y, p, dydt):
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
    __insulin(t, y, p, dydt)
    __fattyacids(t, y, p, dydt)
    __aminoacids(t, y, p, dydt)
    __g6p(t, y, p, dydt)
    __triglycerides(t, y, p, dydt)
    __pyruvate(t, y, p, dydt)
    __acetylcoa(t, y, p, dydt)
    __ROS(t, y, p, dydt)

def __glucose(t, y, p, dydt):
    dydt[Index.plasma_glucose] += (
        + (-p.Subq.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.Subq.k_G_to_plasma * y[Index.subq_glucose] * p.V.subq) / p.V.plasma
        + (-p.Vsc.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.Vsc.k_G_to_plasma * y[Index.vsc_glucose] * p.V.vsc) / p.V.plasma
    )
    dydt[Index.subq_glucose] += (
        + (p.Subq.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.Subq.k_G_to_plasma * y[Index.subq_glucose] * p.V.subq) / p.V.subq
        - p.Subq.k_G_to_G6P * y[Index.subq_glucose] + p.Subq.k_G6P_to_G * y[Index.subq_G6P]
    )
    dydt[Index.vsc_glucose] += (
        + (p.Vsc.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.Vsc.k_G_to_plasma * y[Index.vsc_glucose] * p.V.vsc) / p.V.vsc
        - p.Vsc.k_G_to_G6P * y[Index.vsc_glucose] + p.Vsc.k_G6P_to_G * y[Index.vsc_G6P]
    )


def __insulin(t, y, p, dydt):
    dydt[Index.plasma_insulin] += (
        + (
            - p.Subq.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma
            + p.Subq.k_insulin_to_plasma * y[Index.subq_insulin] * p.V.subq
        ) / p.V.plasma
        + (
            - p.Vsc.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma
            + p.Vsc.k_insulin_to_plasma * y[Index.vsc_insulin] * p.V.vsc
        ) / p.V.plasma
    )
    
    dydt[Index.subq_insulin] += (
        + (p.Subq.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma - p.Subq.k_insulin_to_plasma * y[Index.subq_insulin] * p.V.subq) / p.V.subq
        - p.Subq.kCL_insulin * y[Index.subq_insulin]
    )
    dydt[Index.vsc_insulin] += (
        + (p.Vsc.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma - p.Vsc.k_insulin_to_plasma * y[Index.vsc_insulin] * p.V.vsc) / p.V.vsc
        - p.Vsc.kCL_insulin * y[Index.vsc_insulin]
    )

def __fattyacids(t, y, p, dydt):
    Km = 1
    
    dydt[Index.plasma_fattyacid] += (
        + (
            - p.Subq.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma
            + p.Subq.k_FA_to_plasma * y[Index.subq_fattyacid] * p.V.subq
        ) / p.V.plasma
        + (
            - p.Vsc.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma
            + p.Vsc.k_FA_to_plasma * y[Index.vsc_fattyacid] * p.V.vsc
        ) / p.V.plasma
    )

    dydt[Index.subq_fattyacid] += (
        + (p.Subq.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.Subq.k_FA_to_plasma * y[Index.subq_fattyacid] * p.V.subq) / p.V.subq
        - p.Subq.k_FA_to_ACoA * y[Index.subq_fattyacid]
        - 3 * (p.Subq.k_FA_to_TAG * p.V.subq * y[Index.subq_fattyacid] / (Km + y[Index.subq_fattyacid] * p.V.subq))**3
        + 3 * p.Subq.k_TAG_to_FA * y[Index.subq_TAG]
        + p.Subq.k_ACoA_to_FA * y[Index.subq_ACoA]
    )

    dydt[Index.vsc_fattyacid] += (
        + (p.Vsc.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.Vsc.k_FA_to_plasma * y[Index.vsc_fattyacid] * p.V.vsc) / p.V.vsc
        - p.Subq.k_FA_to_ACoA * y[Index.vsc_fattyacid]
        - 3 * (p.Vsc.k_FA_to_TAG * p.V.vsc * y[Index.vsc_fattyacid] / (Km + y[Index.vsc_fattyacid] * p.V.vsc))**3
        + 3 * p.Vsc.k_TAG_to_FA * y[Index.vsc_TAG]
        + p.Subq.k_ACoA_to_FA * y[Index.vsc_ACoA]
    )



def __aminoacids(t, y, p, dydt):
    dydt[Index.plasma_aminoacid] += (
        + (
            - p.Subq.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            + p.Subq.k_AA_to_plasma * y[Index.subq_aminoacid] * p.V.subq
        ) / p.V.plasma
        + (
            - p.Vsc.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            + p.Vsc.k_AA_to_plasma * y[Index.vsc_aminoacid] * p.V.vsc
        ) / p.V.plasma
    )

    dydt[Index.subq_aminoacid] += (
        + (
            + p.Subq.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            - p.Subq.k_AA_to_plasma * y[Index.subq_aminoacid] * p.V.subq
        ) / p.V.subq
        - p.Subq.k_AA_to_ACoA * y[Index.subq_aminoacid]
    )

    dydt[Index.vsc_aminoacid] += (
        + (
            + p.Vsc.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            - p.Vsc.k_AA_to_plasma * y[Index.vsc_aminoacid] * p.V.vsc
        ) / p.V.vsc
        - p.Subq.k_AA_to_ACoA * y[Index.vsc_aminoacid]
    )


def __g6p(t, y, p, dydt):
    dydt[Index.subq_G6P] += (
        + p.Subq.k_G_to_G6P * y[Index.subq_glucose]
        - p.Subq.k_G6P_to_G * y[Index.subq_G6P]
        - p.Subq.k_G6P_to_P * y[Index.subq_G6P]
        + p.Subq.k_P_to_G6P * y[Index.subq_pyruvate]**2
    )
    dydt[Index.vsc_G6P] += (
        + p.Subq.k_G_to_G6P * y[Index.vsc_glucose]
        - p.Subq.k_G6P_to_G * y[Index.vsc_G6P]
        - p.Subq.k_G6P_to_P * y[Index.vsc_G6P]
        + p.Subq.k_P_to_G6P * y[Index.vsc_pyruvate]**2
    )


def __triglycerides(t, y, p, dydt):
    Km = 1

    dydt[Index.subq_TAG] += (
        + (p.Subq.k_FA_to_TAG * p.V.subq * y[Index.subq_fattyacid] / (Km + y[Index.subq_fattyacid] * p.V.subq))**3
        - p.Subq.k_TAG_to_FA * y[Index.subq_TAG]
    )
    dydt[Index.vsc_TAG] += (
        + (p.Vsc.k_FA_to_TAG * p.V.vsc * y[Index.vsc_fattyacid] / (Km + y[Index.vsc_fattyacid] * p.V.vsc))**3
        - p.Vsc.k_TAG_to_FA * y[Index.vsc_TAG]
    )


def __pyruvate(t, y, p, dydt):
    dydt[Index.subq_pyruvate] += (
        + 2 * p.Subq.k_G6P_to_P * y[Index.subq_G6P]
        - 2 * p.Subq.k_P_to_G6P * y[Index.subq_pyruvate]**2
        - p.Subq.k_P_to_ACoA * y[Index.subq_pyruvate]
    )

    dydt[Index.vsc_pyruvate] += (
        + 2 * p.Subq.k_G6P_to_P * y[Index.vsc_G6P]
        - 2 * p.Subq.k_P_to_G6P * y[Index.vsc_pyruvate]**2
        - p.Subq.k_P_to_ACoA * y[Index.vsc_pyruvate]
    )



def __acetylcoa(t, y, p, dydt):
    dydt[Index.subq_ACoA] += (
        + p.Subq.k_P_to_ACoA * y[Index.subq_pyruvate]
        + 8 * p.Subq.k_FA_to_ACoA * y[Index.subq_fattyacid]
        + p.Subq.k_AA_to_ACoA * y[Index.subq_aminoacid]
        - 8 * p.Subq.k_ACoA_to_FA * y[Index.subq_ACoA]
    )
    dydt[Index.vsc_ACoA] += (
        + p.Subq.k_P_to_ACoA * y[Index.vsc_pyruvate]
        + 8 * p.Subq.k_FA_to_ACoA * y[Index.vsc_fattyacid]
        + p.Subq.k_AA_to_ACoA * y[Index.vsc_aminoacid]
        - 8 * p.Subq.k_ACoA_to_FA * y[Index.vsc_ACoA]
    )

def __ROS(t, y, p, dydt):
    # Fraction of reactions resulting in ROS (1% defa__ult)
    ROSpercent = 0.01
    Km = 1

    dydt[Index.subq_ROS] += ROSpercent * (
        p.Subq.k_FA_to_ACoA * y[Index.subq_fattyacid]
        + p.Subq.k_AA_to_ACoA * y[Index.subq_aminoacid]
        + 3 * (p.Subq.k_FA_to_TAG * p.V.subq * y[Index.subq_fattyacid] / (Km + y[Index.subq_fattyacid] * p.V.subq))**3
        + 3 * p.Subq.k_TAG_to_FA * y[Index.subq_TAG]
        + 8 * p.Subq.k_ACoA_to_FA * y[Index.subq_ACoA]
    )

    dydt[Index.vsc_ROS] += ROSpercent * (
        p.Subq.k_FA_to_ACoA * y[Index.vsc_fattyacid]
        + p.Subq.k_AA_to_ACoA * y[Index.vsc_aminoacid]
        + 3 * (p.Vsc.k_FA_to_TAG * p.V.vsc * y[Index.vsc_fattyacid] / (Km + y[Index.vsc_fattyacid] * p.V.vsc))**3
        + 3 * p.Vsc.k_TAG_to_FA * y[Index.vsc_TAG]
        + 8 * p.Subq.k_ACoA_to_FA * y[Index.vsc_ACoA]
    )
