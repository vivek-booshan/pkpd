from dataclasses import dataclass

# TODO (VIVEK): Actually implement equations in subq, vsc, and muscle that change volume
@dataclass(slots=True)
class Volumes:
    """
    Data class representing the volumes (in liters) of various anatomical compartments.

    Attributes:
        plasma (float): Volume of blood plasma.
        gut (float): Volume of the gastrointestinal tract compartment.
        liver (float): Volume of the liver compartment.
        subq (float): Volume of subcutaneous adipose tissue.
        vsc (float): Volume of visceral adipose tissue.
        muscle (float): Volume of skeletal muscle tissue.
        pancreas (float): Volume of the pancreas.
        brain (float): Volume of the brain.

    Used for scaling concentrations to amounts and for calculating transport dynamics 
    between compartments in the model.
    """
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
    """
    Data class representing shared metabolic rate constants across different tissues.

    All rates are unitless or assumed to be in consistent time units (e.g., 1/min).

    Attributes:
        k_P_to_ACoA (float): Rate constant for conversion of pyruvate (P) to acetyl-CoA (ACoA).
        k_ACoA_to_P (float): Rate constant for conversion of acetyl-CoA (ACoA) back to pyruvate (P).
        k_FA_to_ACoA (float): Rate constant for conversion of fatty acids (FA) to acetyl-CoA (ACoA).
        k_ACoA_to_FA (float): Rate constant for conversion of acetyl-CoA (ACoA) to fatty acids (FA).
        k_AA_to_ACoA (float): Rate constant for conversion of amino acids (AA) to acetyl-CoA (ACoA).
        k_G_to_G6P (float): Rate constant for conversion of glucose (G) to glucose-6-phosphate (G6P).
        k_G6P_to_G (float): Rate constant for conversion of glucose-6-phosphate (G6P) back to glucose (G).
        k_P_to_G6P (float): Rate constant for conversion of pyruvate (P) to glucose-6-phosphate (G6P).
        k_G6P_to_P (float): Rate constant for conversion of glucose-6-phosphate (G6P) to pyruvate (P).

    These parameters are shared across multiple organ models to ensure consistent 
    representation of core metabolic pathways.
    """
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
    """
    Data class representing muscle-specific metabolic rate constants and parameters.

    These parameters govern the dynamic exchanges of various molecules between the muscle and plasma, 
    as well as the muscle's energy metabolism processes.

    Attributes:
        k_insulin_from_plasma (float): Rate constant for the transport of insulin from plasma to muscle.
        k_insulin_to_plasma (float): Rate constant for the transport of insulin from muscle to plasma.
        k_FA_from_plasma (float): Rate constant for the transport of fatty acids from plasma to muscle.
        k_FA_to_plasma (float): Rate constant for the transport of fatty acids from muscle to plasma.
        k_G_from_plasma (float): Rate constant for the transport of glucose from plasma to muscle.
        k_G_to_plasma (float): Rate constant for the transport of glucose from muscle to plasma.
        k_AA_from_plasma (float): Rate constant for the transport of amino acids from plasma to muscle.
        k_AA_to_plasma (float): Rate constant for the transport of amino acids from muscle to plasma.
        k_L_from_plasma (float): Rate constant for the transport of lactate from plasma to muscle.
        k_L_to_plasma (float): Rate constant for the transport of lactate from muscle to plasma.
        NADH_ETC (float): Rate constant for NADH utilization in the electron transport chain (ETC).
        FADH2_ETC (float): Rate constant for FADH2 utilization in the electron transport chain (ETC).
        k_Glc_to_G6P (float): Rate constant for conversion of glucose (Glc) to glucose-6-phosphate (G6P) in muscle.
        k_G6P_to_Glc (float): Rate constant for conversion of glucose-6-phosphate (G6P) to glucose (Glc) in muscle.
        k_P_to_L (float): Rate constant for conversion of pyruvate (P) to lactate (L) in muscle.
        k_L_to_P (float): Rate constant for conversion of lactate (L) to pyruvate (P) in muscle.
        kCL_insulin (float): Clearance rate constant for insulin in muscle.
        kCL_ATP (float): Clearance rate constant for ATP in muscle.
        kCL_FA (float): Clearance rate constant for fatty acids in muscle.

    These parameters are specific to muscle metabolism, regulating nutrient uptake and utilization, 
    and contribute to understanding muscle energy dynamics.
    """
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
    k_ACoA_to_TCA        : float

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
class FatParameters:
    """
    Data class representing parameters governing fatty acid and insulin metabolism in adipose tissue.

    These parameters control the transport and conversion of metabolites between plasma and adipose tissue,
    as well as the regulation of fatty acid storage.

    Attributes:
        k_insulin_from_plasma (float): Rate constant for the transport of insulin from plasma to adipose tissue.
        k_insulin_to_plasma (float): Rate constant for the transport of insulin from adipose tissue to plasma.
        k_FA_from_plasma (float): Rate constant for the transport of fatty acids from plasma to adipose tissue.
        k_FA_to_plasma (float): Rate constant for the transport of fatty acids from adipose tissue to plasma.
        k_G_from_plasma (float): Rate constant for the transport of glucose from plasma to adipose tissue.
        k_G_to_plasma (float): Rate constant for the transport of glucose from adipose tissue to plasma.
        k_AA_from_plasma (float): Rate constant for the transport of amino acids from plasma to adipose tissue.
        k_AA_to_plasma (float): Rate constant for the transport of amino acids from adipose tissue to plasma.
        k_FA_to_TAG (float): Rate constant for conversion of fatty acids to triacylglycerols (TAG) in adipose tissue.
        k_TAG_to_FA (float): Rate constant for conversion of triacylglycerols (TAG) to fatty acids (FA) in adipose tissue.
        kCL_insulin (float): Clearance rate constant for insulin in adipose tissue.

    These parameters regulate the storage and release of fatty acids and other key metabolites in adipose tissue,
    critical for energy storage and utilization.
    """
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
class GIParameters:
    """
    Data class representing gastrointestinal (GI) transport and absorption parameters.

    These parameters govern the absorption, transport, and clearance of nutrients in the gastrointestinal system.

    Attributes:
        kabs_glucose (float): Rate constant for glucose absorption in the gastrointestinal tract.
        kabs_fructose (float): Rate constant for fructose absorption in the gastrointestinal tract.
        kabs_fattyacid (float): Rate constant for fatty acid absorption in the gastrointestinal tract.
        k_diffusion_micelle_to_membrane (float): Rate constant for fatty acid diffusion from micelle to cell membrane.
        k_Vmax_trans (float): Maximum rate of transport of nutrients across the intestinal membrane.
        k_Vmax_reester (float): Maximum rate of re-esterification of fatty acids in the intestinal cells.
        k_Vmax_export (float): Maximum rate of export of absorbed nutrients from the intestinal cells.
        Km_trans (float): Michaelis-Menten constant for nutrient transport across the intestinal membrane.
        Km_reester (float): Michaelis-Menten constant for fatty acid re-esterification in the intestinal cells.
        Km_export (float): Michaelis-Menten constant for nutrient export from the intestinal cells.
        kCL_glucose (float): Clearance rate constant for glucose in the gastrointestinal system.
        kCL_fructose (float): Clearance rate constant for fructose in the gastrointestinal system.
        kCL_fattyacid (float): Clearance rate constant for fatty acids in the gastrointestinal system.

    These parameters are critical for understanding nutrient absorption dynamics, transport efficiency, and clearance 
    from the gastrointestinal tract.
    """
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
class LiverParameters:
    k_FA_from_plasma      : float
    k_FA_to_plasma        : float
    k_G_from_plasma       : float
    k_G_to_plasma         : float
    k_F_from_plasma       : float
    k_F_to_plasma         : float
    k_AA_from_plasma      : float
    k_AA_to_plasma        : float
    k_FA_to_TAG           : float
    k_TAG_to_FA           : float

    k_P_to_ACoA : float
    k_ACoA_to_P : float
    k_FA_to_ACoA: float
    k_ACoA_to_FA: float
    k_AA_to_ACoA: float
    k_G_to_G6P  : float
    k_G6P_to_G  : float
    k_P_to_G6P  : float
    k_G6P_to_P  : float
    pass

@dataclass(frozen=True, slots=True, kw_only=True)
class Parameters:
    """
    Data class representing the parameters used in the metabolic model.

    This class groups together various parameters that control the dynamics of the system, 
    including volumes of different compartments, rate constants for metabolic processes, 
    and the specific parameters governing muscle, adipose tissue (fat), and gastrointestinal (GI) dynamics.

    Attributes:
        V (Volumes): A `Volumes` object containing the volumes of different biological compartments 
                      (e.g., plasma, gut, liver, subq, muscle, etc.).
        Shared (SharedRates): A `SharedRates` object containing rate constants that govern common 
                              metabolic reactions across compartments (e.g., acetyl-CoA conversion rates).
        M (MuscleParameters): A `MuscleParameters` object containing rate constants specific to 
                               muscle metabolism (e.g., glucose, fatty acid, and amino acid transport).
        Subq (FatParameters): A `FatParameters` object containing rate constants specific to 
                              subcutaneous fat metabolism (e.g., fatty acid uptake and conversion rates).
        Vsc (FatParameters): A `FatParameters` object containing rate constants specific to 
                              visceral fat metabolism, similar to the subcutaneous fat parameters.
        GI (GIParameters): A `GIParameters` object containing rate constants governing nutrient absorption 
                            and transport in the gastrointestinal system (e.g., glucose, fructose, and fatty acid transport).

    This class serves as the central container for all parameter values used in the metabolic model, 
    enabling the easy passing and organization of these values to the model's functions and solvers.
    """
    V      : Volumes
    # Shared : SharedRates
    M      : MuscleParameters
    Subq   : FatParameters
    Vsc    : FatParameters
    Liver  : LiverParameters
    GI     : GIParameters

if __name__ == "__main__":
    pass
