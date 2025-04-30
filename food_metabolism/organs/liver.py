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

    __glycogen_synthesis(t, y, p, dydt)
    __glycogen_breakdown(t, y, p, dydt)

    __fattyacid_synthesis(t, y, p, dydt)
    __fattyacid_breakdown(t, y, p, dydt)

    # __aminoacid_metabolism(t, y, p, dydt)

    __insulin(t, y, p, dydt)
    __glucagon(t, y, p, dydt)    
    __somatostatin(t, y, p, dydt)

    # __ATP_synthesis(t, y, p, dydt)
    # 
    __lactate_metabolism(t, y, p, dydt) # plasma lactate uptake and metabolism

    # __cholesterol_synthesis(t, y, p, dydt)
    # __cholesterol_breakdown(t, y, p, dydt)
    #
    # TODO (Vivek): the other infinite number of processes in the liver
    return

def __glucose_uptake(t, y, p, dydt):
    dydt[Index.plasma_glucose] = - y[Index.plasma_glucose] * p.V.plasma * p.Liver.k_G_from_plasma / p.V.liver
    dydt[Index.liver_glucose] = + y[Index.plasma_glucose] * p.V.plasma * p.Liver.k_G_from_plasma / p.V.liver
    return

def __glucose_export(t, y, p, dydt):
    dydt[Index.liver_glucose] = - y[Index.liver_glucose] * p.V.liver * p.Liver.k_G_to_plasma / p.V.plasma
    dydt[Index.plasma_glucose] = + y[Index.liver_glucose] * p.V.liver * p.Liver.k_G_to_plasma / p.V.plasma
    return

def __glycolysis(t, y, p, dydt):
    # NOTE (Vivek): not actually a legit thing, this needs to be refined with better parameters
    # this is just a temporary placeholder to give an idea of what the hell is happening
    dydt[Index.liver_glucose] = - y[Index.liver_glucose] * p.Liver.k_glycolysis
    dydt[Index.liver_extracellular_pyruvate] = + 2 * y[Index.liver_glucose] * p.Liver.k_glycolysis
    dydt[Index.liver_ATP] = + 2 * y[Index.liver_glucose] * p.Liver.k_glycolysis
    dydt[Index.liver_NADH] = + 2 * y[Index.liver_glucose] * p.Liver.k_glycolysis
    dydt[Index.liver_NAD] = - 2 * y[Index.liver_glucose] * p.Liver.k_glycolysis
    return

def __pyruvate_to_mitochondria(t, y, p, dydt):
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

def __lactate_metabolism(t, y, p, dydt):
    __cori_cycle(t, y, p, dydt)
    __lactate_to_pyruvate(t, y, p, dydt)
    return

def __fatty_acid_synthesis(t, y, p, dydt):
    pass
def __fatty_acid_breakdown(t, y, p, dydt):
    pass
