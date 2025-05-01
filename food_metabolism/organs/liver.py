# TODO (Vivek): due to the large number of functions used, cython migration is
# absolutely pertinent to remove the function call costs via static inlining
import numpy as np
from organs.parameters import Parameters
from organs.index import Index

def liver(t: float, y: np.ndarray, p: Parameters) -> np.ndarray:
    dydt = np.zeros(len(Index))
    __liver(t, y, dydt)
    return dydt

def __liver(t: float, y: np.ndarray, p: Parameters, dydt: np.ndarray):
    __glucose_uptake(t, y, p, dydt)
    __glucose_export(t, y, p, dydt)
    __glycolysis(t, y, p, dydt)
    __gluconeogenesis(t, y, p, dydt)

    # __fructose_uptake(t, y, p, dydt)
    # __fructose_export(t, y, p, dydt)
    # __fructose_to_pyruvate(t, y, p, dydt)
    # __pyruvate_to_fructose(t, y, p, dydt)

    # __glycogen_synthesis(t, y, p, dydt)
    # __glycogen_breakdown(t, y, p, dydt)

    # __fattyacid_uptake(t, y, p, dydt)
    # __fattyacid_export(t, y, p, dydt)

    # __aminoacid_uptake(t, y, p, dydt)
    # __aminoacid_export(t, y, p, dydt)

    __insulin_breakdown(t, y, p, dydt)
    __glucagon_breakdown(t, y, p, dydt)    
    __somatostatin_breakdown(t, y, p, dydt)

    __transport_pyruvate_to_mitochondria(t, y, p, dydt)
    __mitochondrial_pyruvate_to_ACoA(t, y, p, dydt)
    __TCA(t, y, p, dydt) 

    # __cholesterol_synthesis(t, y, p, dydt)
    # __cholesterol_breakdown(t, y, p, dydt)
    #
    # TODO (Vivek): the other infinite number of processes in the liver
    return

def __glucose_uptake(t, y, p, dydt):
    dydt[Index.plasma_glucose] += - y[Index.plasma_glucose] * p.V.plasma * p.Liver.k_G_from_plasma / p.V.liver
    dydt[Index.liver_glucose] += + y[Index.plasma_glucose] * p.V.plasma * p.Liver.k_G_from_plasma / p.V.liver
    return

def __glucose_export(t, y, p, dydt):
    dydt[Index.liver_glucose] += - y[Index.liver_glucose] * p.V.liver * p.Liver.k_G_to_plasma / p.V.plasma
    dydt[Index.plasma_glucose] += + y[Index.liver_glucose] * p.V.liver * p.Liver.k_G_to_plasma / p.V.plasma
    return

def __glucose_to_g6p(t, y, p, dydt):
    dydt[Index.liver_G6P] += p.Liver.k_G_to_G6P * y[Index.liver_glucose] + p.Liver.k_G6P_to_G * y[Index.liver_G6P]

    dydt[Index.liver_glucose] += - p.Liver.k_G_to_G6P * y[Index.liver_glucose] + p.Liver.k_G6P_to_G * y[Index.liver_G6P]
    dydt[Index.liver_ATP] += - p.Liver.k_G_to_G6P * y[Index.plasma_insulin] * y[Index.liver_ATP]
    return

def __g6p_to_glucose(t, y, p, dydt):
    dydt[Index.liver_G6P] += - p.Liver.k_G6P_to_G * y[Index.liver_G6P]

    dydt[Index.liver_glucose] += + p.Liver.k_G6P_to_G * y[Index.liver_G6P]
    dydt[Index.liver_ATP] += + p.Liver.k_G6P_to_G * y[Index.liver_G6P]
    return

# NOTE : need to double check this
def __fructose_to_pyruvate(t, y, p, dydt):
    dydt[Index.liver_fructose] += - p.Liver.k_F_to_P * y[Index.liver_fructose] * y[Index.liver_NAD]**2
    dydt[Index.liver_NAD] += - 2 * p.Liver.k_F_to_P * y[Index.liver_G6P] * y[Index.liver_NAD]**2

    dydt[Index.liver_extracellular_pyruvate] += + 2 * p.Liver.k_F_to_P * y[Index.liver_G6P] * y[Index.liver_NAD]**2
    dydt[Index.liver_NADH] += + 2 * p.Liver.k_F_to_P * y[Index.liver_G6P] * y[Index.liver_NAD]**2
    dydt[Index.liver_ATP] += + 3 * p.Liver.k_F_to_P * y[Index.liver_G6P] * y[Index.liver_NAD]**2
    return

# NOTE : need to doublecheck this
def __pyruvate_to_fructose(t, y, p, dydt):
    return

def __pyruvate_to_g6p(t, y, p, dydt):
    dydt[Index.liver_G6P] += + p.Liver.k_P_to_G6P * y[Index.liver_extracellular_pyruvate]**2 * y[Index.liver_ATP]**3 * y[Index.liver_NADH]**2
    dydt[Index.liver_NAD] += + 2 * p.Liver.k_P_to_G6P * y[Index.liver_extracellular_pyruvate]**2 * y[Index.liver_ATP]**3 * y[Index.liver_NADH]**2

    dydt[Index.liver_extracellular_pyruvate] += - 2 * p.Liver.k_P_to_G6P  * y[Index.liver_extracellular_pyruvate]**2 * y[Index.liver_ATP]**3 * y[Index.liver_NADH]**2
    dydt[Index.liver_NADH] += - 2 * p.Liver.k_P_to_G6P * y[Index.liver_extracellular_pyruvate]**2 * y[Index.liver_ATP]**3 * y[Index.liver_NADH]**2
    dydt[Index.liver_ATP] += - 3 * p.Liver.k_P_to_G6P * y[Index.liver_extracellular_pyruvate]**2 * y[Index.liver_ATP]**3 * y[Index.liver_NADH]**2
    return

def __glycolysis(t, y, p, dydt):
    __glucose_to_g6p(t, y, p, dydt)
    __pyruvate_to_g6p(t, y, p, dydt)
    return

def __gluconeogenesis(t, y, p, dydt):
    __pyruvate_to_g6p(t, y, p, dydt)
    __g6p_to_glucose(t, y, p, dydt)
    return

def __fructose_uptake(t, y, p, dydt):
    dydt[Index.plasma_fructose] = - y[Index.plasma_fructose] * p.V.plasma * p.Liver.k_F_from_plasma / p.V.liver
    dydt[Index.liver_fructose] = + y[Index.plasma_fructose] * p.V.plasma * p.Liver.k_F_from_plasma / p.V.liver
    return

def __fructose_export(t, y, p, dydt):
    dydt[Index.liver_fructose] = - y[Index.plasma_fructose] * p.V.plasma * p.Liver.k_F_to_plasma / p.V.liver
    dydt[Index.plasma_fructose] = + y[Index.plasma_fructose] * p.V.plasma * p.Liver.k_F_to_plasma / p.V.liver
    return

def __fructolysis(t, y, p, dydt):
    __fructose_to_pyruvate(t, y, p, dydt)
    return

def __fructoneogenesis(t, y, p, dydt):
    __pyruvate_to_fructose(t, y, p, dydt)
    return

def __transport_pyruvate_to_mitochondria(t, y, p, dydt):
    dydt[Index.liver_extracellular_pyruvate] = - y[Index.liver_extracellular_pyruvate] * p.Liver.k_pyruvate_to_mitochondria
    dydt[Index.liver_mitochondrial_pyruvate] = + y[Index.liver_extracellular_pyruvate] * p.Liver.k_pyruvate_to_mitochondria
    return

def __mitochondrial_pyruvate_to_ACoA(t, y, p, dydt):
    dydt[Index.liver_mitochondrial_pyruvate] = -  p.Liver.k_mitochondrial_pyruvate_to_ACoA * y[Index.liver_mitochondrial_pyruvate] * y[Index.liver_NAD]
    dydt[Index.liver_mitochondrial_ACoA] = + p.Liver.k_mitochondrial_pyruvate_to_ACoA * y[Index.liver_mitochondrial_pyruvate] * y[Index.liver_NAD] 
    dydt[Index.liver_NAD] = - p.Liver.k_mitochondrial_pyruvate_to_ACoA * y[Index.liver_mitochondrial_pyruvate] * y[Index.liver_NAD]
    dydt[Index.liver_NADH] = + p.Liver.k_mitochondrial_pyruvate_to_ACoA * y[Index.liver_mitochondrial_pyruvate] * y[Index.liver_NAD]
    return

def __TCA(t, y, p, dydt):
    # + 2 atp/gtp, 6 NADH, 2 FADH2 per 2 ACoA
    dydt[Index.liver_mitochondrial_ACoA] = - p.Liver.k_ACoA_to_TCA * y[Index.liver_mitochondrial_ACoA] * y[Index.liver_NAD]**3 * y[Index.liver_FAD]
    dydt[Index.liver_NAD] = - 3 * p.Liver.k_ACoA_to_TCA * y[Index.liver_mitochondrial_ACoA] * y[Index.liver_NAD]**3 * y[Index.liver_FAD]
    dydt[Index.liver_NADH] = + 3 * p.Liver.k_ACoA_to_TCA * y[Index.liver_mitochondrial_ACoA] * y[Index.liver_NAD]**3 * y[Index.liver_FAD]
    dydt[Index.liver_FAD] = - p.Liver.k_ACoA_to_TCA * y[Index.liver_mitochondrial_ACoA] * y[Index.liver_NAD]**3 * y[Index.liver_FAD]
    dydt[Index.liver_FADH2] = + p.Liver.k_ACoA_to_TCA * y[Index.liver_mitochondrial_ACoA] * y[Index.liver_NAD]**3 * y[Index.liver_FAD]

    dydt[Index.liver_ATP] = + p.Liver.k_ACoA_to_TCA * y[Index.liver_mitochondrial_ACoA] * y[Index.liver_NAD]**3 * y[Index.liver_FAD]
    return

def __insulin_breakdown(t, y, p, dydt):
    # assume any insulin into liver auto metabolized thus bypass insulin from plasma to liver ecs
    dydt[Index.plasma_insulin] = - p.Liver.kCL_insulin * dydt[Index.plasma_insulin]
    return

def __insulin_modulate(t, y, p, dydt):
    # vary metabolic parameters based on insulin concentration
    pass
   
def __glucagon_breakdown(t, y, p, dydt):
    # assume any glucagon into liver auto metabolized thus bypass glucagon from plasma to liver ecs
    dydt[Index.plasma_glucagon] = - p.Liver.kCL_glucagon * dydt[Index.plasma_glucagon]
    return

def __glucagon_modulate(t, y, p, dydt):
    # vary metabolic parameters based on glucagon concentration
    pass
    
def __somatostatin_breakdown(t, y, p, dydt):
    # assume any somatostatin into liver auto metabolized thus bypass somatostatin from plasma to liver ecs
    dydt[Index.plasma_somatostatin] = - p.Liver.kCL_somatostatin * dydt[Index.plasma_somatostatin]
    return

def __somatostatin_modulate(t, y, p, dydt):
    # vary metabolic parameters based on somatostatin concentration
    pass

def __lactate_uptake(t, y, p, dydt):
    dydt[Index.plasma_lactate] += - p.Liver.k_L_from_plasma * y[Index.plasma_lactate] * p.V.plasma / p.V.liver
    dydt[Index.liver_lactate] += + p.Liver.k_L_from_plasma * y[Index.plasma_lactate] * p.V.plasma / p.V.liver
    return

def __lactate_to_pyruvate(t, y, p, dydt):
    dydt[Index.liver_lactate] += - p.Liver.k_L_to_P * y[Index.liver_lactate] * y[Index.liver_NAD]
    dydt[Index.liver_extracellular_pyruvate] += + p.Liver.k_L_to_P * y[Index.liver_lactate] * y[Index.liver_NAD]
    dydt[Index.liver_NAD] += + p.Liver.k_L_to_P * y[Index.liver_lactate] * y[Index.liver_NAD]
    return

def __ROS(t, y, p, dydt):
    dydt[Index.liver_ROS] += ROSpercent * (
        p.Subq.k_FA_to_ACoA * y[Index.liver_fattyacid]
        + p.Subq.k_AA_to_ACoA * y[Index.liver_aminoacid]
        + 3 * (p.liver.k_FA_to_TAG * p.V.liver * y[Index.liver_fattyacid] / (Km + y[Index.liver_fattyacid] * p.V.liver))**3
        + 3 * p.liver.k_TAG_to_FA * y[Index.liver_TAG]
        + 8 * p.Subq.k_ACoA_to_FA * y[Index.liver_mitochondrial_ACoA]
    )

def __fattyacid_uptake(t, y, p, dydt):
    dydt[Index.plasma_fattyacid] += - p.Liver.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma / p.V.liver
    dydt[Index.liver_fattyacid] += + p.Liver.k_FA_from_plasma * y[Index.plasma_fattyacid] * p.V.plasma / p.V.liver
    return

def __fattyacid_export(t, y, p, dydt):
    dydt[Index.liver_fattyacid] += - p.Liver.k_FA_to_plasma * y[Index.liver_fattyacid] * p.V.liver / p.V.plasma
    dydt[Index.plasma_fattyacid] += + p.Liver.k_FA_to_plasma * y[Index.liver_fattyacid] * p.V.liver / p.V.plasma
    return

def __fattyacid_synthesis(t, y, p, dydt):
    dydt[Index.liver_fattyacid] += - p.Liver.k_FA_to_ACoA * y[Index.liver_fattyacid] * y[Index.liver_ATP]
    dydt[Index.liver_mitochondrial_ACoA] += + 8 * p.Liver.k_FA_to_ACoA * y[Index.liver_fattyacid] * y[Index.liver_ATP]
    return

def __aminoacid_uptake(t, y, p, dydt):
    dydt[Index.plasma_aminoacid] += - p.Liver.k_AA_from_plasma * y[Index.plasma_aminoacid] * p.V.plasma / p.V.liver
    dydt[Index.liver_aminoacid] += + p.Liver.k_AA_to_plasma * y[Index.plasma_aminoacid] * p.V.plasma / p.V.liver
    return

def __aminoacid_export(t, y, p, dydt):
    dydt[Index.plasma_aminoacid] += + p.Liver.k_AA_to_plasma * y[Index.liver_aminoacid] * p.V.liver / p.V.plasma
    dydt[Index.liver_aminoacid] += - p.Liver.k_AA_to_plasma * y[Index.liver_aminoacid] * p.V.liver / p.V.plasma
    return

def __fatty_acid_breakdown(t, y, p, dydt):
    return
