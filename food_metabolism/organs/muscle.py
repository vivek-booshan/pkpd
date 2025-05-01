# TODO (Vivek) : Refactor this file in terms of physiological functions instead of components
# ie using __glycogen_synthesis instead of __glucose and __glycogen
import numpy as np
from .parameters import *
from .index import Index

def skeletalmuscle(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
    """
    The skeletalmuscle function computes the rate of change (dydt) for various metabolites 
    and processes in the skeletal muscle compartment of a metabolic model.

    This function takes in the current state vector, `y`, representing concentrations or amounts 
    of metabolites, and updates it with the rate of change (`dydt`) based on metabolic interactions 
    governed by the parameters `p`.

    Processes:
        - Models glucose, fatty acid, and amino acid uptake, utilization, and storage in skeletal muscle.
        - Simulates the conversion of glucose to glucose-6-phosphate (G6P) and the formation of glycogen.
        - Models the glycolysis and oxidative phosphorylation pathways with the generation of ATP, NADH, FADH2, and ROS.
        - Computes the production of pyruvate and its conversion to acetyl-CoA.
        - Simulates the dynamic formation and utilization of lactate in the muscle.

    Attributes:
        t (float): The current time point (in the simulation's time units).
        y (np.ndarray): The state vector representing the current concentrations or amounts of metabolites in skeletal muscle.
        p (object): A parameters object containing rate constants and volume information.
        n (int): The total number of state variables in the system (e.g., metabolites, molecules involved in muscle metabolism).


    Returns:
        dydt (np.ndarray): The array of rate-of-change values for each state variable, which is used in 
                            numerical integration methods (e.g., `solve_ivp`) to simulate the system.
    """
    dydt = np.zeros(len(Index))
    __muscle(t, y, p, dydt)   
    return dydt

def __muscle(t, y, p, dydt):
    """
    The skeletalmuscle function computes the rate of change (dydt) for various metabolites 
    and processes in the skeletal muscle compartment of a metabolic model.

    This function takes in the current state vector, `y`, representing concentrations or amounts 
    of metabolites, and updates it with the rate of change (`dydt`) based on metabolic interactions 
    governed by the parameters `p`.

    Processes:
        - Models glucose, fatty acid, and amino acid uptake, utilization, and storage in skeletal muscle.
        - Simulates the conversion of glucose to glucose-6-phosphate (G6P) and the formation of glycogen.
        - Models the glycolysis and oxidative phosphorylation pathways with the generation of ATP, NADH, FADH2, and ROS.
        - Computes the production of pyruvate and its conversion to acetyl-CoA.
        - Simulates the dynamic formation and utilization of lactate in the muscle.

    Attributes:
        t (float): The current time point (in the simulation's time units).
        y (np.ndarray): The state vector representing the current concentrations or amounts of metabolites in skeletal muscle.
        p (object): A parameters object containing rate constants and volume information.
        n (int): The total number of state variables in the system (e.g., metabolites, molecules involved in muscle metabolism).
        dydt (np.ndarray): The array of rate-of-change values for each state variable, which is used in 
                            numerical integration methods (e.g., `solve_ivp`) to simulate the system.
    """
    __glucose(t,y,p, dydt)
    __insulin(t,y,p, dydt)
    __fattyacids(t,y,p, dydt)
    __aminoacids(t,y,p, dydt)
    __g6p(t,y,p, dydt)
    __glycogen(t,y,p, dydt)
    __pyruvate(t,y,p, dydt)
    __acetylcoa(t,y,p, dydt)
    __NAD(t,y,p, dydt)
    __NADH(t,y,p, dydt)
    __FAD(t,y,p, dydt)
    __FADH2(t,y,p, dydt)
    __ROS(t,y,p, dydt)
    __ATP(t,y,p, dydt)
    __lactate(t,y,p, dydt)

def __glucose(t, y, p, dydt):
    dydt[Index.plasma_glucose] += (
        (-p.M.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.M.k_G_to_plasma * y[Index.muscle_glucose] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_glucose] += (
        (p.M.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.M.k_G_to_plasma * y[Index.muscle_glucose] * p.V.muscle) / p.V.muscle
        - p.M.k_Glc_to_G6P * y[Index.muscle_glucose] + p.M.k_G6P_to_Glc * y[Index.muscle_G6P]
    )

def __insulin(t, y, p, dydt):
    dydt[Index.plasma_insulin] += (
        (-p.M.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma + p.M.k_insulin_to_plasma * y[Index.muscle_insulin] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_insulin] += (
        (p.M.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma - p.M.k_insulin_to_plasma * y[Index.muscle_insulin] * p.V.muscle) / p.V.muscle
        - p.M.kCL_insulin * y[Index.muscle_insulin]
    )

def __fattyacids(t, y, p, dydt):
    dydt[Index.plasma_fattyacid] += (
        (-p.M.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma + p.M.k_FA_to_plasma * y[Index.muscle_fattyacid] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_fattyacid] += (
        (p.M.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.M.k_FA_to_plasma * y[Index.muscle_fattyacid] * p.V.muscle) / p.V.muscle
        - p.M.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
    )

def __aminoacids(t, y, p, dydt):
    dydt[Index.plasma_aminoacid] += (
        (-p.M.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma + p.M.k_AA_to_plasma * y[Index.muscle_aminoacid] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_aminoacid] += (
        (p.M.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma - p.M.k_AA_to_plasma * y[Index.muscle_aminoacid] * p.V.muscle) / p.V.muscle
        - p.M.k_AA_to_ACoA * y[Index.muscle_aminoacid]
    )

def __g6p(t, y, p, dydt):
    Km = 1
    dydt[Index.muscle_G6P] += (
        + p.M.k_Glc_to_G6P * y[Index.muscle_glucose]
        - p.M.k_G6P_to_Glc * y[Index.muscle_G6P]
        - p.M.k_G6P_to_P * p.V.muscle * y[Index.muscle_G6P] / (Km + y[Index.muscle_G6P] * p.V.muscle)
        + p.M.k_P_to_G6P * y[Index.muscle_glycogen]
        - p.M.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        + p.M.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
    )

def __glycogen(t, y, p, dydt):
    #adult weighing 70kg has about 400 g glycogenin muscle
    #1-2% muscle mass
    #20% of volume
    Km = 1

    dydt[Index.muscle_glycogen] += (
        p.M.k_G6P_to_P * p.V.muscle * y[Index.muscle_G6P] / (Km + y[Index.muscle_G6P] * p.V.muscle)
        - p.M.k_P_to_G6P * y[Index.muscle_glycogen]
    )

def __pyruvate(t, y, p, dydt):
    dydt[Index.muscle_pyruvate] += (
        2 * p.M.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        - 2 * p.M.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        - p.M.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        - p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH]
        + p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
    )

def __acetylcoa(t, y, p, dydt):
    dydt[Index.muscle_ACoA] += (
        p.M.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        + 8 * p.M.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
        + p.M.k_AA_to_ACoA * y[Index.muscle_aminoacid]
        - p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        - p.M.k_ACoA_to_TCA * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD] 
    )

def __NAD(t, y, p, dydt):
    dydt[Index.muscle_NAD] += (
        -2 * p.M.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        + 2 * p.M.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        - p.M.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        - 3 * p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + 2 * p.M.NADH_ETC * y[Index.muscle_NADH]**2
        + p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH]
        - p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
        - 3 * p.M.k_ACoA_to_TCA * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD] 
    )

def __NADH(t, y, p, dydt):
    dydt[Index.muscle_NADH] += (
        2 * p.M.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        - 2 * p.M.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        + p.M.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        + 3 * p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        - 2 * p.M.NADH_ETC * y[Index.muscle_NADH]**2
        - p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH]
        + p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
        + 3 * p.M.k_ACoA_to_TCA * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD] 
    )

def __FAD(t, y, p, dydt):
    dydt[Index.muscle_FAD] += (
        -p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + 2 * p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
        - p.M.k_ACoA_to_TCA * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD] 
    )

def __FADH2(t, y, p, dydt):
    dydt[Index.muscle_FADH2] += (
        p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        - 2 * p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
        + p.M.k_ACoA_to_TCA * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD] 
    )

def __ROS(t, y, p, dydt):
    # might include a randomness aspect here. not sure what the probability will be
    ROSpercent = 0.02
    # 1-2% of molecular oxygen is converted to superoxide owing to electron leak
    dydt[Index.muscle_ROS] += ROSpercent * (
        p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
        + p.M.NADH_ETC * y[Index.muscle_NADH]**2
        + p.M.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
        + p.M.k_AA_to_ACoA * y[Index.muscle_aminoacid]
    )


def __ATP(t, y, p, dydt):
    dydt[Index.muscle_ATP] += (
        -p.M.k_G_to_G6P * y[Index.plasma_insulin] * y[Index.muscle_ATP]
        + p.M.k_G6P_to_G * y[Index.muscle_G6P]
        + 3 * p.M.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        - 3 * p.M.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        + p.M.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + 3 * p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
        + 5 * p.M.NADH_ETC * y[Index.muscle_NADH]**2
        - p.M.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
        - p.M.kCL_ATP * y[Index.muscle_ATP]
    )

def __lactate(t, y, p, dydt):
    dydt[Index.plasma_lactate] += (
        (-p.M.k_L_from_plasma * y[Index.plasma_lactate] * p.V.plasma + p.M.k_L_to_plasma * y[Index.muscle_lactate] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_lactate] += (
        (p.M.k_L_from_plasma * y[Index.plasma_lactate] * p.V.plasma - p.M.k_L_to_plasma * y[Index.muscle_lactate] * p.V.muscle) / p.V.muscle
        + p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH] - p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
    )


