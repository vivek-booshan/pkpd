from dataclasses import dataclass
from enum import IntEnum, auto

class Index(IntEnum):
    plasma_glucose=auto()
    plasma_fructose=auto()
    plasma_insulin=auto()
    plasma_fattyacid=auto()
    plasma_aminoacid=auto()
    plasma_lactate=auto()

    subq_glucose=auto()
    subq_insulin=auto()
    subq_fattyacid=auto()
    subq_aminoacid=auto()
    subq_G6P=auto()
    subq_TAG=auto()
    subq_pyruvate=auto()
    subq_ACoA=auto()
    subq_ROS=auto()

    vsc_glucose=auto()
    vsc_insulin=auto()
    vsc_fattyacid=auto()
    vsc_aminoacid=auto()
    vsc_G6P=auto()
    vsc_TAG=auto()
    vsc_pyruvate=auto()
    vsc_ACoA=auto()
    vsc_ROS=auto()

    muscle_glucose=auto()
    muscle_insulin=auto()
    muscle_fattyacid=auto()
    muscle_aminoacid=auto()
    muscle_G6P=auto()
    muscle_glycogen=auto()
    muscle_pyruvate=auto()
    muscle_ACoA=auto()
    muscle_NAD=auto()
    muscle_NADH=auto()
    muscle_FAD=auto()
    muscle_FADH2=auto()
    muscle_ROS=auto()
    muscle_ATP=auto()
    muscle_lactate=auto()

    gut_glucose=auto()
    gut_fructose=auto()
    micellar_fattyacid=auto()
    membrane_fattyacid=auto()
    cytosol_fattyacid=auto()
    cytosol_TAG=auto()

def get_plasma_indices(index: Index):
    return tuple(
        index.plasma_glucose,
        index.plasma_insulin,
        index.plasma_fattyacid,
        index.plasma_aminoacid,
        index.plasma_lactate,
    )

def get_subq_indices(index: Index):
    return tuple(
        index.subq_glucose,
        index.subq_insulin,
        index.subq_fattyacid,
        index.subq_aminoacid,
        index.subq_G6P,
        index.subq_TAG,
        index.subq_ACoA,
        index.subq_pyruvate,
        index.subq_ROS,
    )

def get_vsc_indices(index: Index):
    return tuple(
        index.vsc_glucose,
        index.vsc_insulin,
        index.vsc_fattyacid,
        index.vsc_aminoacid,
        index.vsc_G6P,
        index.vsc_TAG,
        index.vsc_ACoA,
        index.vsc_pyruvate,
        index.vsc_ROS,
    )

def get_muscle_indices(index: Index):
    return tuple(
        index.muscle_glucose,
        index.muscle_insulin,
        index.muscle_fattyacid,
        index.muscle_aminoacid,
        index.muscle_G6P,
        index.muscle_glycogen,
        index.muscle_pyruvate,
        index.muscle_AcoA,
        index.muscle_NAD,
        index.muscle_NADH,
        index.muscle_FAD,
        index.muscle_FADH2,
        index.muscle_ROS,
        index.muscle_ATP,
        index.muscle_lactate,
    )

def get_gut_indices(index: Index):
    return tuple(
        index.gut_glucose,
        index.gut_fructose,
        index.micellar_fattyacid,
        index.membrane_fattyacid,
        index.cytosol_fattyacid,
        index.cytosol_TAG,
    )



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
class PancreasParameters:
    pass

@dataclass(frozen=True, slots=True, kw_only=True)
class Parameters:
    V      : Volumes
    shared : SharedRates
    CL     : ClearanceRates
    M      : MuscleParameters
    Subq   : FatParameters
    Vsc    : FatParameters
    GI     : GIParameters

if __name__ == "__main__":
    pass
