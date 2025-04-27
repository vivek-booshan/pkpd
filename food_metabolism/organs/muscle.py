import numpy as np
from organs.parameters import *
from organs.index import Index

def skeletalmuscle(t, y, p):
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

    glucose(t,y,p, dydt)
    insulin(t,y,p, dydt)
    fattyacids(t,y,p, dydt)
    aminoacids(t,y,p, dydt)
    g6p(t,y,p, dydt)
    glycogen(t,y,p, dydt)
    pyruvate(t,y,p, dydt)
    acetylcoa(t,y,p, dydt)
    NAD(t,y,p, dydt)
    NADH(t,y,p, dydt)
    FAD(t,y,p, dydt)
    FADH2(t,y,p, dydt)
    ROS(t,y,p, dydt)
    ATP(t,y,p, dydt)
    lactate(t,y,p, dydt)
    
    return dydt

def __muscle(t, y, p, dydt):
    glucose(t,y,p, dydt)
    insulin(t,y,p, dydt)
    fattyacids(t,y,p, dydt)
    aminoacids(t,y,p, dydt)
    g6p(t,y,p, dydt)
    glycogen(t,y,p, dydt)
    pyruvate(t,y,p, dydt)
    acetylcoa(t,y,p, dydt)
    NAD(t,y,p, dydt)
    NADH(t,y,p, dydt)
    FAD(t,y,p, dydt)
    FADH2(t,y,p, dydt)
    ROS(t,y,p, dydt)
    ATP(t,y,p, dydt)
    lactate(t,y,p, dydt)

def glucose(t, y, p, dydt):
    dydt[Index.plasma_glucose] += (
        (-p.M.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.M.k_G_to_plasma * y[Index.muscle_glucose] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_glucose] += (
        (p.M.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.M.k_G_to_plasma * y[Index.muscle_glucose] * p.V.muscle) / p.V.muscle
        - p.M.k_Glc_to_G6P * y[Index.muscle_glucose] + p.M.k_G6P_to_Glc * y[Index.muscle_G6P]
    )

def insulin(t, y, p, dydt):
    dydt[Index.plasma_insulin] += (
        (-p.M.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma + p.M.k_insulin_to_plasma * y[Index.muscle_insulin] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_insulin] += (
        (p.M.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma - p.M.k_insulin_to_plasma * y[Index.muscle_insulin] * p.V.muscle) / p.V.muscle
        - p.M.kCL_insulin * y[Index.muscle_insulin]
    )

def fattyacids(t, y, p, dydt):
    dydt[Index.plasma_fattyacid] += (
        (-p.M.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma + p.M.k_FA_to_plasma * y[Index.muscle_fattyacid] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_fattyacid] += (
        (p.M.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.M.k_FA_to_plasma * y[Index.muscle_fattyacid] * p.V.muscle) / p.V.muscle
        - p.Shared.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
    )

def aminoacids(t, y, p, dydt):
    dydt[Index.plasma_aminoacid] += (
        (-p.M.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma + p.M.k_AA_to_plasma * y[Index.muscle_aminoacid] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_aminoacid] += (
        (p.M.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma - p.M.k_AA_to_plasma * y[Index.muscle_aminoacid] * p.V.muscle) / p.V.muscle
        - p.Shared.k_AA_to_ACoA * y[Index.muscle_aminoacid]
    )

def g6p(t, y, p, dydt):
    Km = 1
    dydt[Index.muscle_G6P] += (
        + p.M.k_Glc_to_G6P * y[Index.muscle_glucose]
        - p.M.k_G6P_to_Glc * y[Index.muscle_G6P]
        - p.Shared.k_G6P_to_P * p.V.muscle * y[Index.muscle_G6P] / (Km + y[Index.muscle_G6P] * p.V.muscle)
        + p.Shared.k_P_to_G6P * y[Index.muscle_glycogen]
        - p.Shared.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        + p.Shared.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
    )

def glycogen(t, y, p, dydt):
    #adult weighing 70kg has about 400 g glycogenin muscle
    #1-2% muscle mass
    #20% of volume
    Km = 1

    dydt[Index.muscle_glycogen] += (
        p.Shared.k_G6P_to_P * p.V.muscle * y[Index.muscle_G6P] / (Km + y[Index.muscle_G6P] * p.V.muscle)
        - p.Shared.k_P_to_G6P * y[Index.muscle_glycogen]
    )

def pyruvate(t, y, p, dydt):
    dydt[Index.muscle_pyruvate] += (
        2 * p.Shared.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        - 2 * p.Shared.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        - p.Shared.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        - p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH]
        + p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
    )

def acetylcoa(t, y, p, dydt):
    dydt[Index.muscle_ACoA] += (
        p.Shared.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        + 8 * p.Shared.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
        + p.Shared.k_AA_to_ACoA * y[Index.muscle_aminoacid]
        - p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
    )

def NAD(t, y, p, dydt):
    dydt[Index.muscle_NAD] += (
        -2 * p.Shared.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        + 2 * p.Shared.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        - p.Shared.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        - 3 * p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + 2 * p.M.NADH_ETC * y[Index.muscle_NADH]**2
        + p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH]
        - p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
    )

def NADH(t, y, p, dydt):
    dydt[Index.muscle_NADH] += (
        2 * p.Shared.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        - 2 * p.Shared.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        + p.Shared.k_P_to_ACoA * y[Index.muscle_pyruvate] * y[Index.muscle_NAD]
        + 3 * p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        - 2 * p.M.NADH_ETC * y[Index.muscle_NADH]**2
        - p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH]
        + p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
    )

def FAD(t, y, p, dydt):
    dydt[Index.muscle_FAD] += (
        -p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + 2 * p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
    )

def FADH2(t, y, p, dydt):
    dydt[Index.muscle_FADH2] += (
        p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        - 2 * p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
    )

def ROS(t, y, p, dydt):
    # might include a randomness aspect here. not sure what the probability will be
    ROSpercent = 0.02
    # 1-2% of molecular oxygen is converted to superoxide owing to electron leak
    dydt[Index.muscle_ROS] += ROSpercent * (
        p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
        + p.M.NADH_ETC * y[Index.muscle_NADH]**2
        + p.Shared.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
        + p.Shared.k_AA_to_ACoA * y[Index.muscle_aminoacid]
    )


def ATP(t, y, p, dydt):
    dydt[Index.muscle_ATP] += (
        -p.Shared.k_G_to_G6P * y[Index.plasma_insulin] * y[Index.muscle_ATP]
        + p.Shared.k_G6P_to_G * y[Index.muscle_G6P]
        + 3 * p.Shared.k_G6P_to_P * y[Index.muscle_G6P] * y[Index.muscle_NAD]**2
        - 3 * p.Shared.k_P_to_G6P * y[Index.muscle_pyruvate]**2 * y[Index.muscle_ATP]**3 * y[Index.muscle_NADH]**2
        + p.Shared.k_ACoA_to_P * y[Index.muscle_ACoA] * y[Index.muscle_NAD]**3 * y[Index.muscle_FAD]
        + 3 * p.M.FADH2_ETC * y[Index.muscle_FADH2]**2
        + 5 * p.M.NADH_ETC * y[Index.muscle_NADH]**2
        - p.Shared.k_FA_to_ACoA * y[Index.muscle_fattyacid] * y[Index.muscle_ATP]
        - p.M.kCL_ATP * y[Index.muscle_ATP]
    )

def lactate(t, y, p, dydt):
    dydt[Index.plasma_lactate] += (
        (-p.M.kCL_FA * y[Index.plasma_lactate] * p.V.plasma + p.M.kCL_FA * y[Index.muscle_lactate] * p.V.muscle) / p.V.plasma
    )
    dydt[Index.muscle_lactate] += (
        (p.M.kCL_FA * y[Index.plasma_lactate] * p.V.plasma - p.M.kCL_FA * y[Index.muscle_lactate] * p.V.muscle) / p.V.muscle
        + p.M.k_P_to_L * y[Index.muscle_pyruvate] * y[Index.muscle_NADH] - p.M.k_L_to_P * y[Index.muscle_lactate] * y[Index.muscle_NAD]
    )

def muscle_init():
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            gut=0.0,
            liver=0.0,
            subq=0.0,
            vsc=0.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
        ),
        Shared=SharedRates(
            k_P_to_ACoA=1,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=5,               # acetylcoa_to_TCA
            k_FA_to_ACoA=1/8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1/4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=0.0,            # not provided
            k_G_to_G6P=1,                # glucose_to_g6p
            k_G6P_to_G=0.1,              # g6p_to_glucose
            k_P_to_G6P=0.1,              # pyruvate_to_g6p
            k_G6P_to_P=1                 # g6p_to_pyruvate
        ),
        M=MuscleParameters(
            k_insulin_from_plasma=5,
            k_insulin_to_plasma=0.5,
            k_FA_from_plasma=1,
            k_FA_to_plasma=0.1,
            k_G_from_plasma=1,
            k_G_to_plasma=0.1,
            k_AA_from_plasma=1,
            k_AA_to_plasma=0.1,
            k_L_from_plasma=0.1,         # lactate_plasma_skeletalmuscle
            k_L_to_plasma=1,             # lactate_skeletalmuscle_plasma
            NADH_ETC=1,
            FADH2_ETC=1,
            k_Glc_to_G6P=1,              # glucose_to_g6p (duplicate of Shared, but muscle-local here)
            k_G6P_to_Glc=0.1,            # g6p_to_glucose
            k_P_to_L=0.1,                # pyruvate_to_lactate
            k_L_to_P=0.1,                # lactate_to_pyruvate
            kCL_insulin=1,
            kCL_ATP=1,
            kCL_FA=1,
        ),
        Subq=None,
        Vsc=None,
        GI=None
    )
    return p
