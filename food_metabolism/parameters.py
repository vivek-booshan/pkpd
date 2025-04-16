from dataclasses import dataclass

@dataclass(slots=True)
class Volumes:
    plasma  : float
    gut     : float
    liver   : float
    fat     : float
    muscle  : float
    pancreas: float
    brain   : float

@dataclass(frozen=True, slots=True)
class SharedRates:
    k_P_to_ACoA : float
    k_ACoA_to_P : float
    k_FA_to_ACoA: float
    k_AA_to_ACoA: float
    k_AcoA_to_FA: float
    k_G_to_G6P  : float
    k_G6P_to_G  : float
    k_P_to_G6P  : float
    k_G6P_to_P  : float

@dataclass(frozen=True, slots=True)
class ClearanceRates:
    kCL_insulin: float
    kCL_ATP    : float
    kCL_G      : float
    kCL_F      : float
    kCL_FA     : float

@dataclass(frozen=True, slots=True)
class MuscleParameters:
    k_insulin_from_plasma: float
    k_insulin_to_plasma  : float
    k_FA_from_plasma     : float
    k_FA_to_plasma       : float
    k_G_from_plasma      : float
    k_G_to_plasma        : float
    k_AA_from_plasma     : float
    k_AA_to_plasma       : float
    k_L_from_plasma      : float
    k_L_to_plasma        : float
    NADH_ETC             : float
    FADH2_ETC            : float
    k_Glc_to_G6P         : float
    k_G6P_to_Glc         : float
    k_P_to_L             : float
    k_L_to_P             : float

@dataclass(frozen=True, slots=True)
class FatParameters:
    k_insulin_from_plasma : float
    k_insulin_to_plasma   : float
    k_FA_from_plasma      : float
    k_FA_to_plasma        : float
    k_G_from_plasma       : float
    k_G_to_plasma         : float
    k_AA_from_plasma      : float
    k_AA_to_plasma        : float
    k_FA_to_TAG           : float
    k_TAG_to_FA           : float

@dataclass(frozen=True, slots=True)
class GIParameters:
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

@dataclass(frozen=True, slots=True)
class PancreasParameters(frozen=True, slots=True):
    pass

@dataclass(slots=True)
class Parameters:
    V      : Volumes
    shared : SharedRates
    CL     : ClearanceRates
    M      : MuscleParameters
    F      : FatParameters
    GI     : GIParameters
