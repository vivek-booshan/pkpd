from dataclasses import dataclass

@dataclass(slots=True)
class Volumes:
    plasma  : float
    gut     : float
    liver   : float
    subq    : float
    vsc     : float
    muscle  : float
    pancreas: float
    brain   : float

@dataclass(frozen=True, slots=True)
class SharedRates:
    k_P_to_ACoA : float
    k_ACoA_to_P : float
    k_FA_to_ACoA: float
    k_ACoA_to_FA: float
    k_AA_to_ACoA: float
    k_G_to_G6P  : float
    k_G6P_to_G  : float
    k_P_to_G6P  : float
    k_G6P_to_P  : float

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
    kCL_insulin          : float
    kCL_ATP              : float
    kCL_FA               : float

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
    kCL_insulin           : float

@dataclass(frozen=True, slots=True)
class GIParameters:
    kabs_glucose                   : float
    kabs_fructose                  : float
    kabs_fattyacid                 : float
    k_diffusion_micelle_to_membrane: float
    k_Vmax_trans                   : float
    k_Vmax_reester                 : float
    k_Vmax_export                  : float
    Km_trans                       : float
    Km_reester                     : float
    Km_export                      : float
    kCL_glucose                    : float
    kCL_fructose                   : float
    kCL_fattyacid                  : float

@dataclass(frozen=True, slots=True)
class PancreasParameters:
    pass

@dataclass(frozen=True, slots=True, kw_only=True)
class Parameters:
    V      : Volumes
    Shared : SharedRates
    M      : MuscleParameters
    Subq   : FatParameters
    Vsc    : FatParameters
    GI     : GIParameters

if __name__ == "__main__":
    pass
