# 0 glucose plasma
# 1 glucose fat
# 2 insulin plasma
# 3 insulin fat
# 4 fattyacids plasma
# 5 fattyacids fat
# 6 aminoacids plasma
# 7 aminoacids fat
# 8 g6p fat
# 9 triglycerides fat
# 10 Pyruvate
# 11 acetylcoa fat
# 12 ROS fat
import numpy as np
from parameters import *

def fat_init():
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            fat=11.0,
            gut=0.0,
            liver=0.0,
            muscle=0.0,
            pancreas=0.0,
            brain=0.0
        ),
        shared=SharedRates(
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
        CL=ClearanceRates(
            kCL_insulin=1.0,
            kCL_ATP=0.0,
            kCL_G=0.0,
            kCL_F=0.0,
            kCL_FA=0.0
        ),
        F=FatParameters(
            k_insulin_from_plasma=1.0,
            k_insulin_to_plasma=0.1,
            k_FA_from_plasma=2.0,
            k_FA_to_plasma=0.2,
            k_G_from_plasma=1.0,
            k_G_to_plasma=0.1,
            k_AA_from_plasma=1.0,
            k_AA_to_plasma=0.1,
            k_FA_to_TAG=1.0,               # fattyacids_to_triglycerides
            k_TAG_to_FA=0.1                # triglycerides_to_fattyacids
        ),
        M=None,
        GI=None,
    )
    return p


def fat(t, y, p, n):
    dydt = np.zeros(n)
    glucose(t, y, p, dydt)
    insulin(t, y, p, dydt)
    fattyacids(t, y, p, dydt)
    aminoacids(t, y, p, dydt)
    g6p(t, y, p, dydt)
    triglycerides(t, y, p, dydt)
    pyruvate(t, y, p, dydt)
    acetylcoa(t, y, p, dydt)
    ROS(t, y, p, dydt)

    # dydt = dglucose + dinsulin + dfattyacids + daminoacids + dg6p + dtriglycerides + dpyruvate + dacetylcoa + dROS
    return dydt


def glucose(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[0] = (
        (-p.F.k_G_from_plasma * y[0] * p.V.plasma + p.F.k_G_to_plasma * y[1] * p.V.fat) / p.V.plasma
    )
    dydt[1] = (
        (p.F.k_G_from_plasma * y[0] * p.V.plasma - p.F.k_G_to_plasma * y[1] * p.V.fat) / p.V.fat
        - p.shared.k_G_to_G6P * y[1] + p.shared.k_G6P_to_G * y[8]
    )
    # return dydt


def insulin(t, y, p, dydt):
    # dydt = np.zeros(n)
    dydt[2] = (
        (-p.F.k_insulin_from_plasma * y[2] * p.V.plasma + p.F.k_insulin_to_plasma * y[3] * p.V.fat) / p.V.plasma
    )
    dydt[3] = (
        + (p.F.k_insulin_from_plasma * y[2] * p.V.plasma - p.F.k_insulin_to_plasma * y[3] * p.V.fat) / p.V.fat
        - p.CL.kCL_insulin * y[3]
    )
    # return dydt

def fattyacids(t, y, p, dydt):
    Km = 1
    # dydt = np.zeros(n)
    
    dydt[4] = (
        (-p.F.k_FA_from_plasma * y[4] * p.V.plasma + p.F.k_FA_to_plasma * y[5] * p.V.fat) / p.V.plasma
    )
    dydt[5] = (
        + (p.F.k_FA_from_plasma * y[4] * p.V.plasma - p.F.k_FA_to_plasma * y[5] * p.V.fat) / p.V.fat
        - p.shared.k_FA_to_ACoA * y[5]
        - 3 * (p.F.k_FA_to_TAG * p.V.fat * y[5] / (Km + y[5] * p.V.fat))**3
        + 3 * p.F.k_TAG_to_FA * y[9]
        + p.shared.k_ACoA_to_FA * y[11]
    )
    # return dydt


def aminoacids(t, y, p, dydt):
    # dydt = np.zeros(n)

    dydt[6] = (
        (-p.F.k_AA_from_plasma * y[6] * p.V.plasma + p.F.k_AA_to_plasma * y[7] * p.V.fat) / p.V.plasma
    )
    dydt[7] = (
        + (p.F.k_AA_from_plasma * y[6] * p.V.plasma - p.F.k_AA_to_plasma * y[7] * p.V.fat) / p.V.fat
        - p.shared.k_AA_to_ACoA * y[7]
    )
    # return dydt


def g6p(t, y, p, dydt):
    # dydt = np.zeros(n)
    
    dydt[8] = (
        + p.shared.k_G_to_G6P * y[1]
        - p.shared.k_G6P_to_G * y[8]
        - p.shared.k_G6P_to_P * y[8]
        + p.shared.k_P_to_G6P * y[10]**2
    )
    # return dydt


def triglycerides(t, y, p, dydt):
    Km = 1
    # dydt = np.zeros(n)

    dydt[9] = (
        + (p.F.k_FA_to_TAG * p.V.fat * y[5] / (Km + y[5] * p.V.fat))**3
        - p.F.k_TAG_to_FA * y[9]
    )
    # return dydt


def pyruvate(t, y, p, dydt):
    # dydt = np.zeros(n)

    dydt[10] = (
        + 2 * p.shared.k_G6P_to_P * y[8]
        - 2 * p.shared.k_P_to_G6P * y[10]**2
        - p.shared.k_P_to_ACoA * y[10]
    )
    # return dydt


def acetylcoa(t, y, p, dydt):
    # dydt = np.zeros(n)

    dydt[11] = (
        + p.shared.k_P_to_ACoA * y[10]
        + 8 * p.shared.k_FA_to_ACoA * y[5]
        + p.shared.k_AA_to_ACoA * y[7]
        - 8 * p.shared.k_ACoA_to_FA * y[11]
    )
    # return dydt

def ROS(t, y, p, dydt):
    # Fraction of reactions resulting in ROS (1% default)
    ROSpercent = 0.01
    Km = 1

    # dydt = np.zeros(n)

    dydt[12] = ROSpercent * (
        p.shared.k_FA_to_ACoA * y[5]
        + p.shared.k_AA_to_ACoA * y[7]
        + 3 * (p.F.k_FA_to_TAG * p.V.fat * y[5] / (Km + y[5] * p.V.fat))**3
        + 3 * p.F.k_TAG_to_FA * y[9]
        + 8 * p.shared.k_ACoA_to_FA * y[11]
    )

    # return dydt
