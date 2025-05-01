from .fat import __fat
from .gi import __GI
from .muscle import __muscle
from .pancreas import __pancreas
from .liver import __liver
from .index import Index
from .parameters import *

import numpy as np
def system(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
    dydt = np.zeros(len(Index))
    __fat(t, y, p, dydt)
    __GI(t, y, p, dydt)
    __muscle(t, y, p, dydt)
    __pancreas(t, y, dydt)
    __liver(t, y, p, dydt)
    return dydt

def init() -> Parameters:
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            subq=11.0,
            vsc=1.0,
            gut=1.25,
            liver=2.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
        ),
        M=MuscleParameters(
            k_insulin_from_plasma=5,
            k_insulin_to_plasma=0.5,
            k_FA_from_plasma=1,
            k_FA_to_plasma=0.1,
            k_G_from_plasma=1,
            k_G_to_plasma=0.1,
            k_AA_from_plasma=1,
            k_ACoA_to_TCA=10,
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
            kCL_FA=0,
            k_P_to_ACoA=1.0,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=0.0,               # unused
            k_FA_to_ACoA=1 / 8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1 / 4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=0,                # acetylcoa_to_fattyacids
            k_G_to_G6P=1.0,                # glucose_to_g6p
            k_G6P_to_G=0.1,                # g6p_to_glucose
            k_P_to_G6P=0.1,                # pyruvate_to_g6p
            k_G6P_to_P=1.0,                # g6p_to_pyruvate
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
            k_TAG_to_FA=0.01,              # triglycerides_to_fattyacids
            kCL_insulin=1.0,
            k_P_to_ACoA=1.0,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=0.0,               # unused
            k_FA_to_ACoA=1 / 8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1 / 4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=1.0,              # acetylcoa_to_fattyacids
            k_G_to_G6P=1.0,                # glucose_to_g6p
            k_G6P_to_G=0.1,                # g6p_to_glucose
            k_P_to_G6P=0.1,                # pyruvate_to_g6p
            k_G6P_to_P=1.0,                # g6p_to_pyruvate
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
            k_TAG_to_FA=0.01,             # triglycerides_to_fattyacids
            kCL_insulin=1.0,
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
        ),
        Liver=LiverParameters(
            k_G_from_plasma=1,
            k_G_to_plasma=0.1,
            k_F_from_plasma=1,
            k_F_to_plasma=0.1,
            k_FA_from_plasma=1,
            k_FA_to_plasma=0.1,
            k_AA_from_plasma=1,
            k_AA_to_plasma=0.1,

            k_FA_to_TAG=1.0,               # fattyacids_to_triglycerides
            k_TAG_to_FA=0.01,              # triglycerides_to_fattyacids
            kCL_insulin=1.0,
            kCL_glucagon=1,
            kCL_somatostatin=1,
            k_P_to_ACoA=1.0,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=0.0,               # unused
            k_FA_to_ACoA=1 / 8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1 / 4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=1.0,              # acetylcoa_to_fattyacids
            k_G_to_G6P=1.0,                # glucose_to_g6p
            k_G6P_to_G=0.1,                # g6p_to_glucose
            k_P_to_G6P=0.1,                # pyruvate_to_g6p
            k_G6P_to_P=1.0,                # g6p_to_pyruvate

           

        ),
    )
    return p

def giInit() -> Parameters:
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            subq=11.0,
            vsc=1.0,
            gut=1.25,
            liver=0.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
        ),
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
        ),
        M=None,
        Subq=None,
        Vsc=None,
        Liver=None,
    )
    return p

def muscle_init() -> Parameters:
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            subq=11.0,
            vsc=1.0,
            gut=1.25,
            liver=0.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
        ),
        M=MuscleParameters(
            k_insulin_from_plasma=5,
            k_insulin_to_plasma=0.5,
            k_FA_from_plasma=1,
            k_FA_to_plasma=0.1,
            k_G_from_plasma=1,
            k_G_to_plasma=0.1,
            k_AA_from_plasma=1,
            k_ACoA_to_TCA=10,
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
            kCL_FA=0,
            k_P_to_ACoA=1.0,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=0.0,               # unused
            k_FA_to_ACoA=1 / 8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1 / 4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=0,              # acetylcoa_to_fattyacids
            k_G_to_G6P=1.0,                # glucose_to_g6p
            k_G6P_to_G=0.1,                # g6p_to_glucose
            k_P_to_G6P=0.1,                # pyruvate_to_g6p
            k_G6P_to_P=1.0                 # g6p_to_pyruvate
        ),
        Subq=None,
        Vsc=None,
        Liver=None,
        GI=None,
    )
    return p

def subq_init() -> Parameters:
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            subq=11.0,
            vsc=1.0,
            gut=1.25,
            liver=0.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
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
            k_TAG_to_FA=0.01,              # triglycerides_to_fattyacids
            kCL_insulin=1.0,
            k_P_to_ACoA=1.0,               # pyruvate_to_acetylcoa
            k_ACoA_to_P=0.0,               # unused
            k_FA_to_ACoA=1 / 8,            # fattyacids_to_acetylcoa
            k_AA_to_ACoA=1 / 4,            # aminoacids_to_acetylcoa
            k_ACoA_to_FA=1.0,              # acetylcoa_to_fattyacids
            k_G_to_G6P=1.0,                # glucose_to_g6p
            k_G6P_to_G=0.1,                # g6p_to_glucose
            k_P_to_G6P=0.1,                # pyruvate_to_g6p
            k_G6P_to_P=1.0,                # g6p_to_pyruvate
        ),
        M=None,
        Vsc=None,
        Liver=None,
        GI=None,
    )
    return p

def vsc_init() -> Parameters:
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            subq=11.0,
            vsc=1.0,
            gut=1.25,
            liver=0.0,
            muscle=25.0,
            pancreas=0.0,
            brain=0.0
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
            k_TAG_to_FA=0.01,             # triglycerides_to_fattyacids
            kCL_insulin=1.0,
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
        Liver=None,
        M=None,
        GI=None,
    )
    return p
