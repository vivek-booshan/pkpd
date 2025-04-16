import numpy as np

def glucose_two_compartment(t, y, params):
    # TODO : need to change this with respect to a holistic y
    # (ie: F_gut, F_blood = y[m], y[n])
    G_gut, G_blood = y

    V_gut = params.V.gut
    V_blood = params.V.plasma  # Using plasma as "blood" pool

    kabs = params.GI.kabs_G
    kclear = params.CL.kCL_G

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
    kclear = params.CL.kCL_F

    dF_gut_dt = -(kabs * F_gut * V_gut) / V_gut
    dF_blood_dt = ((kabs * F_gut * V_gut) - (kclear * F_blood * V_blood)) / V_blood

    return np.array([dF_gut_dt, dF_blood_dt])

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
    J_clear = CL.kCL_FA * FA_blood

    dFA_mic_dt = -J_diff
    dFA_mem_dt = J_diff - J_trans
    dFA_cyto_dt = J_trans - J_reester
    dTG_cyto_dt = J_reester - J_export
    dFA_blood_dt = (J_export - J_clear)

    return np.array([dFA_mic_dt, dFA_mem_dt, dFA_cyto_dt, dTG_cyto_dt, dFA_blood_dt])
