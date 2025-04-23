import numpy as np
from parameters import *


def skeletalmuscle(t, y, p, n):
    """
    - 0 glucose plasma
    - 1 glucose muscle
    - 2 insulin plasma
    - 3 insulin muscle
    - 4 fattyacids plasma
    - 5 fattyacids muscle
    - 6 aminoacids plasma
    - 7 aminoacids muscle
    - 8 g6p muscle
    - 9 glycogen muscle
    - 10 Pyruvate
    - 11 acetylcoa muscle
    - 12 NAD 
    - 13 NADH
    - 14 FAD
    - 15 FADH2
    - 16 ROS
    - 17 ATP
    - 18 lactate in plasma
    - 19 lactate in muscle
    """
    dydt = np.zeros(n)

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

def glucose(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[0] = (
        (-p.M.k_G_from_plasma * y[0] * p.V.plasma + p.M.k_G_to_plasma * y[1] * p.V.muscle) / p.V.plasma
    )
    dydt[1] = (
        (p.M.k_G_from_plasma * y[0] * p.V.plasma - p.M.k_G_to_plasma * y[1] * p.V.muscle) / p.V.muscle
        - p.M.k_Glc_to_G6P * y[1] + p.M.k_G6P_to_Glc * y[8]
    )
    # return dydt

def insulin(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[2] = (
        (-p.M.k_insulin_from_plasma * y[2] * p.V.plasma + p.M.k_insulin_to_plasma * y[3] * p.V.muscle) / p.V.plasma
    )
    dydt[3] = (
        (p.M.k_insulin_from_plasma * y[2] * p.V.plasma - p.M.k_insulin_to_plasma * y[3] * p.V.muscle) / p.V.muscle
        - p.CL.kCL_insulin * y[3]
    )
    # return dydt

def fattyacids(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[4] = (
        (-p.M.k_FA_from_plasma * y[4] * p.V.plasma + p.M.k_FA_to_plasma * y[5] * p.V.muscle) / p.V.plasma
    )
    dydt[5] = (
        (p.M.k_FA_from_plasma * y[4] * p.V.plasma - p.M.k_FA_to_plasma * y[5] * p.V.muscle) / p.V.muscle
        - p.shared.k_FA_to_ACoA * y[5] * y[17]
    )
    # return dydt

def aminoacids(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[6] = (
        (-p.M.k_AA_from_plasma * y[6] * p.V.plasma + p.M.k_AA_to_plasma * y[7] * p.V.muscle) / p.V.plasma
    )
    dydt[7] = (
        (p.M.k_AA_from_plasma * y[6] * p.V.plasma - p.M.k_AA_to_plasma * y[7] * p.V.muscle) / p.V.muscle
        - p.shared.k_AA_to_ACoA * y[7]
    )
    # return dydt

def g6p(t, y, p, dydt):
    # dydt = np.zeros(n)
    Km = 1
    dydt[8] = (
        + p.M.k_Glc_to_G6P * y[1]
        - p.M.k_G6P_to_Glc * y[8]
        - p.shared.k_G6P_to_P * p.V.muscle * y[8] / (Km + y[8] * p.V.muscle)
        + p.shared.k_P_to_G6P * y[9]
        - p.shared.k_G6P_to_P * y[8] * y[12]**2
        + p.shared.k_P_to_G6P * y[10]**2 * y[17]**3 * y[13]**2
    )
    # return dydt

def glycogen(t, y, p, dydt):
    # dydt = np.zeros(n)
    #adult weighing 70kg has about 400 g glycogenin muscle
    #1-2% muscle mass
    #20% of volume
    Km = 1

    dydt[9] = (
        p.shared.k_G6P_to_P * p.V.muscle * y[8] / (Km + y[8] * p.V.muscle)
        - p.shared.k_P_to_G6P * y[9]
    )
    # return dydt

def pyruvate(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[10] = (
        2 * p.shared.k_G6P_to_P * y[8] * y[12]**2
        - 2 * p.shared.k_P_to_G6P * y[10]**2 * y[17]**3 * y[13]**2
        - p.shared.k_P_to_ACoA * y[10] * y[12]
        - p.M.k_P_to_L * y[10] * y[13]
        + p.M.k_L_to_P * y[19] * y[12]
    )
    # return dydt

def acetylcoa(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[11] = (
        p.shared.k_P_to_ACoA * y[10] * y[12]
        + 8 * p.shared.k_FA_to_ACoA * y[5] * y[17]
        + p.shared.k_AA_to_ACoA * y[7]
        - p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
    )
    # return dydt

def NAD(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[12] = (
        -2 * p.shared.k_G6P_to_P * y[8] * y[12]**2
        + 2 * p.shared.k_P_to_G6P * y[10]**2 * y[17]**3 * y[13]**2
        - p.shared.k_P_to_ACoA * y[10] * y[12]
        - 3 * p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
        + 2 * p.M.NADH_ETC * y[13]**2
        + p.M.k_P_to_L * y[10] * y[13]
        - p.M.k_L_to_P * y[19] * y[12]
    )
    # return dydt

def NADH(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[13] = (
        2 * p.shared.k_G6P_to_P * y[8] * y[12]**2
        - 2 * p.shared.k_P_to_G6P * y[10]**2 * y[17]**3 * y[13]**2
        + p.shared.k_P_to_ACoA * y[10] * y[12]
        + 3 * p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
        - 2 * p.M.NADH_ETC * y[13]**2
        - p.M.k_P_to_L * y[10] * y[13]
        + p.M.k_L_to_P * y[19] * y[12]
    )
    # return dydt

def FAD(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[14] = (
        -p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
        + 2 * p.M.FADH2_ETC * y[15]**2
    )
    # return dydt

def FADH2(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[15] = (
        p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
        - 2 * p.M.FADH2_ETC * y[15]**2
    )
    # return dydt

def ROS(t, y, p, dydt):
    # might include a randomness aspect here. not sure what the probability will be
    ROSpercent = 0.02
    # 1-2% of molecular oxygen is converted to superoxide owing to electron leak
    # dydt = np.zeros(n)
    dydt[16] = ROSpercent * (
        p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
        + p.M.FADH2_ETC * y[15]**2
        + p.M.NADH_ETC * y[13]**2
        + p.shared.k_FA_to_ACoA * y[5] * y[17]
        + p.shared.k_AA_to_ACoA * y[7]
    )
    # return dydt


def ATP(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[17] = (
        -p.shared.k_G_to_G6P * y[2] * y[17]
        + p.shared.k_G6P_to_G * y[8]
        + 3 * p.shared.k_G6P_to_P * y[8] * y[12]**2
        - 3 * p.shared.k_P_to_G6P * y[10]**2 * y[17]**3 * y[13]**2
        + p.shared.k_ACoA_to_P * y[11] * y[12]**3 * y[14]
        + 3 * p.M.FADH2_ETC * y[15]**2
        + 5 * p.M.NADH_ETC * y[13]**2
        - p.shared.k_FA_to_ACoA * y[5] * y[17]
        - p.CL.kCL_ATP * y[17]
    )
    # return dydt

def lactate(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[18] = (
        (-p.CL.kCL_FA * y[18] * p.V.plasma + p.CL.kCL_FA * y[19] * p.V.muscle) / p.V.plasma
    )
    dydt[19] = (
        (p.CL.kCL_FA * y[18] * p.V.plasma - p.CL.kCL_FA * y[19] * p.V.muscle) / p.V.muscle
        + p.M.k_P_to_L * y[10] * y[13] - p.M.k_L_to_P * y[19] * y[12]
    )
    # return dydt

def muscle_init():
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            gut=0.0,
            liver=0.0,
            fat=0.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
        ),
        shared=SharedRates(
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
        CL=ClearanceRates(
            kCL_insulin=1,               # CLinsulin
            kCL_ATP=1,                   # CLATP
            kCL_G=0.0,
            kCL_F=0.0,
            kCL_FA=0.0
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
            k_Glc_to_G6P=1,              # glucose_to_g6p (duplicate of shared, but muscle-local here)
            k_G6P_to_Glc=0.1,            # g6p_to_glucose
            k_P_to_L=0.1,                # pyruvate_to_lactate
            k_L_to_P=0.1                 # lactate_to_pyruvate
        ),
        F=FatParameters(
            k_insulin_from_plasma=0.0,
            k_insulin_to_plasma=0.0,
            k_FA_from_plasma=0.0,
            k_FA_to_plasma=0.0,
            k_G_from_plasma=0.0,
            k_G_to_plasma=0.0,
            k_AA_from_plasma=0.0,
            k_AA_to_plasma=0.0,
            k_FA_to_TAG=0.0,
            k_TAG_to_FA=0.0
        ),
        GI=GIParameters(
            kabs_G=0.0,
            kabs_F=0.0,
            kabs_FA=0.0,
            k_diffusion_micelle_to_membrane=0.0,
            k_Vmax_trans=0.0,
            k_Vmax_reester=0.0,
            k_Vmax_export=0.0,
            Km_trans=0.0,
            Km_reester=0.0,
            Km_export=0.0
        )
    )
    return p