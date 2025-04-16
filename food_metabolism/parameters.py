from dataclasses import dataclass

@dataclass
class Volumes(slots=True):
    plasma  : float
    gut     : float
    liver   : float
    fat     : float
    muscle  : float
    pancreas: float
    brain   : float

@dataclass
class SharedRates(frozen=True, slots=True):
    k_P_to_ACoA : float
    k_ACoA_to_P : float
    k_FA_to_ACoA: float
    k_AA_to_ACoA: float
    k_AcoA_to_FA: float
    k_G_to_G6P  : float
    k_G6P_to_G  : float
    k_P_to_G6P  : float
    k_G6P_to_P  : float

@dataclass
class ClearanceRates(frozen=True, slots=True):
    kCL_insulin: float
    kCL_ATP    : float
    kCL_G      : float
    kCL_F      : float
    kCL_FA     : float

@dataclass
class MuscleParameters(frozen=True, slots=True):
    k_insulin_plasma: float
    k_insulin_muscle: float
    k_FA_plasma     : float
    k_FA_muscle     : float
    k_G_plasma      : float
    k_G_muscle      : float
    k_AA_plasma     : float
    k_AA_muscle     : float
    k_L_plasma      : float
    k_L_muscle      : float
    NADH_ETC        : float
    FADH2_ETC       : float
    k_Glc_to_G6P    : float
    k_G6P_to_Glc    : float
    k_P_to_L        : float
    k_L_to_P        : float

@dataclass
class FatParameters(frozen=True, slots=True):
    k_insulin_plasma: float
    k_insulin_fat   : float
    k_FA_plasma     : float
    k_FA_fat        : float
    k_G_plasma      : float
    k_G_fat         : float
    k_AA_plasma     : float
    k_AA_fat        : float
    k_FA_to_TAG     : float
    k_TAG_to_FA     : float

@dataclass
class GIParameters(frozen=True, slots=True):
    kabs_G                         : float
    kabs_F                         : float
    kabs_FA                        : float
    k_diffusion_micelle_to_membrane: float
    k_Vmax_trans                   : float
    k_Vmax_reester                 : float
    k_Vmax_export                  : float
    Km_trans                       : float
    Km_reester                     : float
    Km_export                      : float

@dataclass
class PancreasParameters(frozen=True, slots=True):
    pass

@dataclass
class Parameters(slots=True):
    V      : Volumes
    shared : SharedRates
    CL     : ClearanceRates
    M      : MuscleParameters
    F      : FatParameters
    GI     : GIParameters
