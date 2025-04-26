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
            subq=11.0,
            vsc=1,
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
        Subq=FatParameters(
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
        Vsc=FatParameters(
            k_insulin_from_plasma=2.0,
            k_insulin_to_plasma=0.1,
            k_FA_from_plasma=4,
            k_FA_to_plasma=0.2,
            k_G_from_plasma=2.0,
            k_G_to_plasma=0.1,
            k_AA_from_plasma=2.0,
            k_AA_to_plasma=0.1,
            k_FA_to_TAG=1.0,               # fattyacids_to_triglycerides
            k_TAG_to_FA=0.0005             # triglycerides_to_fattyacids
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

    return dydt


def glucose(t, y, p, dydt):
    dydt[Index.plasma_glucose] = (
        + (-p.Subq.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.Subq.k_G_to_plasma * y[Index.subq_glucose] * p.V.subq) / p.V.plasma
        + (-p.Vsc.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma + p.Vsc.k_G_to_plasma * y[Index.vsc_glucose] * p.V.vsc) / p.V.plasma
    )
    dydt[Index.subq_glucose] = (
        + (p.Subq.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.Subq.k_G_to_plasma * y[Index.subq_glucose] * p.V.subq) / p.V.subq
        - p.shared.k_G_to_G6P * y[Index.subq_glucose] + p.shared.k_G6P_to_G * y[Index.subq_G6P]
    )
    dydt[Index.vsc_glucose] = (
        + (p.Vsc.k_G_from_plasma * y[Index.plasma_glucose] * p.V.plasma - p.Vsc.k_G_to_plasma * y[Index.vsc_glucose] * p.V.vsc) / p.V.vsc
        - p.shared.k_G_to_G6P * y[Index.vsc_glucose] + p.shared.k_G6P_to_G * y[Index.vsc_G6P]
    )


def insulin(t, y, p, dydt):
    dydt[Index.plasma_insulin] = (
        + (
            - p.Subq.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma
            + p.Subq.k_insulin_to_plasma * y[Index.subq_insulin] * p.V.subq
        ) / p.V.plasma
        + (
            - p.Vsc.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma
            + p.Vsc.k_insulin_to_plasma * y[Index.vsc_insulin] * p.V.vsc
        ) / p.V.plasma
    )
    
    dydt[Index.subq_insulin] = (
        + (p.Subq.k_insulin_from_plasma * y[Index.plasma_insulin] * p.V.plasma - p.Subq.k_insulin_to_plasma * y[Index.subq_insulin] * p.V.subq) / p.V.subq
        - p.CL.kCL_insulin * y[Index.subq_insulin]
    )

def fattyacids(t, y, p, dydt):
    Km = 1
    
    dydt[Index.plasma_fattyacid] = (
        + (
            - p.Subq.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma
            + p.Subq.k_FA_to_plasma * y[Index.subq_fattyacid] * p.V.subq
        ) / p.V.plasma
        + (
            - p.Vsc.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma
            + p.Vsc.k_FA_to_plasma * y[Index.vsc_fattyacid] * p.V.vsc
        ) / p.V.plasma
    )

    dydt[Index.subq_fattyacid] = (
        + (p.Subq.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.Subq.k_FA_to_plasma * y[Index.subq_fattyacid] * p.V.subq) / p.V.subq
        - p.shared.k_FA_to_ACoA * y[Index.subq_fattyacid]
        - 3 * (p.Subq.k_FA_to_TAG * p.V.subq * y[Index.subq_fattyacid] / (Km + y[Index.subq_fattyacid] * p.V.subq))**3
        + 3 * p.Subq.k_TAG_to_FA * y[Index.subq_TAG]
        + p.shared.k_ACoA_to_FA * y[Index.subq_ACoA]
    )

    dydt[Index.vsc_fattyacid] = (
        + (p.Vsc.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma - p.Vsc.k_FA_to_plasma * y[Index.vsc_fattyacid] * p.V.vsc) / p.V.vsc
        - p.shared.k_FA_to_ACoA * y[Index.vsc_fattyacid]
        - 3 * (p.Vsc.k_FA_to_TAG * p.V.vsc * y[Index.vsc_fattyacid] / (Km + y[Index.vsc_fattyacid] * p.V.vsc))**3
        + 3 * p.Vsc.k_TAG_to_FA * y[Index.vsc_TAG]
        + p.shared.k_ACoA_to_FA * y[Index.vsc_ACoA]
    )



def aminoacids(t, y, p, dydt):
    # dydt = np.zeros(n)

    dydt[Index.plasma_aminoacid] = (
        + (
            - p.Subq.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            + p.Subq.k_AA_to_plasma * y[Index.subq_aminoacid] * p.V.subq
        ) / p.V.plasma
        + (
            - p.Vsc.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            + p.Vsc.k_AA_to_plasma * y[Index.vsc_aminoacid] * p.V.vsc
        ) / p.V.plasma
    )

    dydt[Index.subq_aminoacid] = (
        + (
            + p.Subq.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            - p.Subq.k_AA_to_plasma * y[Index.subq_aminoacid] * p.V.subq
        ) / p.V.subq
        - p.shared.k_AA_to_ACoA * y[Index.subq_aminoacid]
    )

    dydt[Index.vsc_aminoacid] = (
        + (
            + p.Vsc.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma
            - p.Vsc.k_AA_to_plasma * y[Index.vsc_aminoacid] * p.V.vsc
        ) / p.V.vsc
        - p.shared.k_AA_to_ACoA * y[Index.vsc_aminoacid]
    )


def g6p(t, y, p, dydt):
    # dydt = np.zeros(n)
    
    dydt[Index.subq_G6P] = (
        + p.shared.k_G_to_G6P * y[Index.subq_glucose]
        - p.shared.k_G6P_to_G * y[Index.subq_G6P]
        - p.shared.k_G6P_to_P * y[Index.subq_G6P]
        + p.shared.k_P_to_G6P * y[Index.subq_pyruvate]**2
    )
    dydt[Index.vsc_G6P] = (
        + p.shared.k_G_to_G6P * y[Index.vsc_glucose]
        - p.shared.k_G6P_to_G * y[Index.vsc_G6P]
        - p.shared.k_G6P_to_P * y[Index.vsc_G6P]
        + p.shared.k_P_to_G6P * y[Index.vsc_pyruvate]**2
    )
    # return dydt


def triglycerides(t, y, p, dydt):
    Km = 1

    dydt[Index.subq_TAG] = (
        + (p.Subq.k_FA_to_TAG * p.V.subq * y[Index.subq_fattyacid] / (Km + y[Index.subq_fattyacid] * p.V.subq))**3
        - p.Subq.k_TAG_to_FA * y[Index.subq_TAG]
    )
    dydt[Index.vsc_TAG] = (
        + (p.Vsc.k_FA_to_TAG * p.V.vsc * y[Index.vsc_fattyacid] / (Km + y[Index.vsc_fattyacid] * p.V.vsc))**3
        - p.Vsc.k_TAG_to_FA * y[Index.vsc_TAG]
    )


def pyruvate(t, y, p, dydt):
    dydt[Index.subq_pyruvate] = (
        + 2 * p.shared.k_G6P_to_P * y[Index.subq_G6P]
        - 2 * p.shared.k_P_to_G6P * y[Index.subq_pyruvate]**2
        - p.shared.k_P_to_ACoA * y[Index.subq_pyruvate]
    )

    dydt[Index.vsc_pyruvate] = (
        + 2 * p.shared.k_G6P_to_P * y[Index.vsc_G6P]
        - 2 * p.shared.k_P_to_G6P * y[Index.vsc_pyruvate]**2
        - p.shared.k_P_to_ACoA * y[Index.vsc_pyruvate]
    )



def acetylcoa(t, y, p, dydt):

    dydt[Index.subq_ACoA] = (
        + p.shared.k_P_to_ACoA * y[Index.subq_pyruvate]
        + 8 * p.shared.k_FA_to_ACoA * y[Index.subq_fattyacid]
        + p.shared.k_AA_to_ACoA * y[Index.subq_aminoacid]
        - 8 * p.shared.k_ACoA_to_FA * y[Index.subq_ACoA]
    )
    dydt[Index.vsc_ACoA] = (
        + p.shared.k_P_to_ACoA * y[Index.vsc_pyruvate]
        + 8 * p.shared.k_FA_to_ACoA * y[Index.vsc_fattyacid]
        + p.shared.k_AA_to_ACoA * y[Index.vsc_aminoacid]
        - 8 * p.shared.k_ACoA_to_FA * y[Index.vsc_ACoA]
    )

def ROS(t, y, p, dydt):
    # Fraction of reactions resulting in ROS (1% default)
    ROSpercent = 0.01
    Km = 1

    dydt[Index.subq_ROS] = ROSpercent * (
        p.shared.k_FA_to_ACoA * y[Index.subq_fattyacid]
        + p.shared.k_AA_to_ACoA * y[Index.subq_aminoacid]
        + 3 * (p.Subq.k_FA_to_TAG * p.V.subq * y[Index.subq_fattyacid] / (Km + y[Index.subq_fattyacid] * p.V.subq))**3
        + 3 * p.Subq.k_TAG_to_FA * y[Index.subq_TAG]
        + 8 * p.shared.k_ACoA_to_FA * y[Index.subq_ACoA]
    )

    dydt[Index.vsc_ROS] = ROSpercent * (
        p.shared.k_FA_to_ACoA * y[Index.vsc_fattyacid]
        + p.shared.k_AA_to_ACoA * y[Index.vsc_aminoacid]
        + 3 * (p.Vsc.k_FA_to_TAG * p.V.vsc * y[Index.vsc_fattyacid] / (Km + y[Index.vsc_fattyacid] * p.V.vsc))**3
        + 3 * p.Vsc.k_TAG_to_FA * y[Index.vsc_TAG]
        + 8 * p.shared.k_ACoA_to_FA * y[Index.vsc_ACoA]
    )
