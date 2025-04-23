import numpy as np
from parameters import *

def giInit():
    p = Parameters(
        V=Volumes(
            plasma=5.0,
            fat=11.0,
            gut=1.25,
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
            kCL_insulin=0,
            kCL_ATP=0.0,
            kCL_G_GI=0.05,
            kCL_F_GI=0.05,
            kCL_FA_GI=0.1
        ),
        F=None,
        M=None,
        GI=GIParameters(
            kabs_G=0.1,
            kabs_F=0.1,
            kabs_FA=0.05,
            k_diffusion_micelle_to_membrane=0.1,
            k_Vmax_trans=5.0,
            k_Vmax_reester=3.0,
            k_Vmax_export=2.0,
            Km_trans=10.0,
            Km_reester=5.0,
            Km_export=4.0,
        )
        # Pancreas is not being used, so not included
    )
    return p

def glucose_two_compartment(t, y, params):
    # TODO : need to change this with respect to a holistic y
    # (ie: F_gut, F_blood = y[m], y[n])
    G_gut, G_blood = y

    V_gut = params.V.gut
    V_blood = params.V.plasma  # Using plasma as "blood" pool

    kabs = params.GI.kabs_G
    kclear = params.CL.kCL_G_GI

    dG_gut_dt = -(kabs * G_gut * V_gut) / V_gut
    dG_blood_dt = ((kabs * G_gut * V_gut) - (kclear * G_blood * V_blood)) / V_blood

    return np.array([dG_gut_dt, dG_blood_dt])

def fructose_two_compartment(t, y, params):
    # TODO : need to change this with respect to a holistic y
    # (ie: F_gut, F_blood = y[m], y[n])
    F_gut, F_blood = y

    V_gut = params.V.gut
    V_blood = params.V.plasma

    kabs = params.GI.kabs_F
    kclear = params.CL.kCL_F_GI

    dF_gut_dt = -(kabs * F_gut * V_gut) / V_gut
    dF_blood_dt = ((kabs * F_gut * V_gut) - (kclear * F_blood * V_blood)) / V_blood

    return np.array([dF_gut_dt, dF_blood_dt])

def fatty_acid_two_compartment(t, y, params):
    FA_gut, FA_blood = y

    V_gut = params.V.gut
    V_blood = params.V.plasma

    kabs = params.GI.kabs_FA
    kclear = params.CL.kCL_FA_GI

    dFA_gut_dt = -(kabs * FA_gut * V_gut) / V_gut
    dFA_blood_dt = ((kabs * FA_gut * V_gut) - (kclear * FA_blood * V_blood)) / V_blood
    return [dFA_gut_dt, dFA_blood_dt]

def fatty_acid_full_model(t, y, params):
    # TODO : need to change this with respect to a holistic y
    # (ie: F_gut, F_blood = y[m], y[n])
    FA_mic, FA_mem, FA_cyto, TG_cyto, FA_blood = y

    GI = params.GI
    CL = params.CL
    V_blood = params.V.plasma

    # Pull out parameters
    J_diff = GI.k_diffusion_micelle_to_membrane * (FA_mic - FA_mem)
    J_trans = (GI.k_Vmax_trans * FA_mem) / (GI.Km_trans + FA_mem + 1e-6)
    J_reester = (GI.k_Vmax_reester * FA_cyto) / (GI.Km_reester + FA_cyto + 1e-6)
    J_export = (GI.k_Vmax_export * TG_cyto) / (GI.Km_export + TG_cyto + 1e-6)
    J_clear = CL.kCL_FA_GI * FA_blood

    dFA_mic_dt = -J_diff
    dFA_mem_dt = J_diff - J_trans
    dFA_cyto_dt = J_trans - J_reester
    dTG_cyto_dt = J_reester - J_export
    dFA_blood_dt = (J_export - J_clear)

    return np.array([dFA_mic_dt, dFA_mem_dt, dFA_cyto_dt, dTG_cyto_dt, dFA_blood_dt])
