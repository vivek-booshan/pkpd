import numpy as np
import random
from dataclasses import dataclass 
import ligand

INSULIN_ACTIVE = False
NORMAL_BLOOD_SUGAR_LOW : float = 4 # mmol/L subject to change
NORMAL_BLOOD_SUGAR_HIGH : float = 7 # mmol/L

@dataclass(frozen=True, slots=True)
class PancreasParameters:
    # INSULIN
    kinsulin_binding_on          : float
    kinsulin_binding_off         : float
    kinsulin_max_secretion       : float
    kinsulin_tuner               : float
    kinsulin_pancreas_portal     : float
    kinsulin_liver_blood         : float

    # TODO : GLUCAGON
    kglucagon_binding_on         : float
    kglucagon_binding_off        : float
    kglucagon_max_secretion      : float
    kglucagon_tuner              : float
    kglucagon_pancreas_portal    : float
    kglucagon_liver_blood        : float
     
    # TODO : SOMATOSTATIN
    ksomatostatin_binding_on     : float
    ksomatostatin_binding_off    : float
    ksomatostatin_max_secretion  : float
    ksomatostatin_tuner          : float
    # do we care about somatostatin effect on rest of body?
    # mainly just attenuates insulin and glucagon
    ksomatostatin_pancreas_portal: float
    ksomatostatin_liver_blood    : float



def pancreas(t, Ins, RIns, Glc, RGlc, Somat, RSomat, Glu, RGlu, p_pancreas: PancreasParameters):
    """
    Main function to calculate the derivative (dydt) for each component.
    Args:
        t: Timestamp (not used in this case, but passed for compatibility)
        Ins: Free Insulin Blood Conc.
        RIns: Available Insulin Receptor Conc.
        Glc: Free Glucagon Blood Conc.
        RGlc: Available Glucagon Receptor Conc.
        Somat: Free Somatostatin Blood Conc.
        RSomat: Available Somatostatin Receptor Conc.
        Glu: Free Glucose Blood Conc.
        RGlu: Available Glucose Receptor Conc.
        p_pancreas: Parameter dataclass containing all required constants
    Returns:
        dydt: A list containing the derivatives of each component
    """
    dinsulin = insulin(Ins, RIns, p_pancreas)
    dglucagon = glucagon(Glc, RGlc, p_pancreas)
    dsomatostatin = somatostatin(Somat, RSomat, p_pancreas)
    dglucose = glucose(Glu, RGlu, p_pancreas)

    dydt = [
        dinsulin,
        dglucagon,
        dsomatostatin,
        dglucose
    ]

    return dydt

def insulin(I, R, p: PancreasParameters):
    """
    Function to calculate the insulin derivative.
    Args:
        I: Free Insulin Blood concentration
        R: Free Insulin Receptor concentration
        p: Parameter dataclass containing constants for the system
    Returns:
        dinsulin: The calculated derivative for insulin
    """
    # Using a placeholder for ligand's methods
    dinsulin = ligand.ligandODE(I, R, p.kinsulin_binding_on, p.kinsulin_binding_off)
    
    # Modify insulin secretion based on inhibitory secretion
    # NOTE : necessary ? 
    dinsulin += ligand.inhibitory_secretion(I, R, p.kinsulin_binding_on, p.kinsulin_binding_off, p.kinsulin_baseline_secretion_inhibitory, p.kinsulin_tuner)
    
    # Other logic to modulate insulin (e.g., upregulation or inhibition by other factors) can be added here
    # dinsulin += glucose_factor  # Example for glucose effect
    # dinsulin -= somatostatin_factor  # Example for somatostatin effect

    return dinsulin

def glucagon(G, R, p: PancreasParameters):
    """
    Function to calculate the glucagon derivative.
    Args:
        G: Free Glucagon Blood concentration
        R: Free Glucagon Receptor concentration
        p: Parameter dataclass containing constants for the system
    Returns:
        dglucagon: The calculated derivative for glucagon
    """
    # Placeholder for glucagon dynamics, assuming some simple equations.
    dglucagon = (p.glucagon_blood_pancreas * R * p.V_blood - p.glucagon_pancreas_portal * R * p.V_pancreas) / p.V_pancreas
    return dglucagon

def somatostatin(Somat, RSomat, p: PancreasParameters):
    """
    Function to calculate the somatostatin derivative.
    Args:
        Somat: Free Somatostatin Blood concentration
        RSomat: Available Somatostatin Receptor concentration
        p: Parameter dataclass containing constants for the system
    Returns:
        dsomatostatin: The calculated derivative for somatostatin
    """
    dsomatostatin = (p.somatostatin_blood_pancreas * RSomat * p.V_blood - p.somatostatin_pancreas_portal * RSomat * p.V_pancreas) / p.V_pancreas
    return dsomatostatin

def glucose(Glu, RGlu, p: PancreasParameters):
    """
    Function to calculate the glucose derivative.
    Args:
        Glu: Free Glucose Blood concentration
        RGlu: Available Glucose Receptor concentration
        p: Parameter dataclass containing constants for the system
    Returns:
        dglucose: The calculated derivative for glucose
    """
    dglucose = (p.glucose_blood_pancreas * RGlu * p.V_blood - p.glucose_pancreas_portal * RGlu * p.V_pancreas) / p.V_pancreas
    return dglucose

def beta_cell(t, y, p):
    dydt = np.zeros((n, 1))
    # Gplasma = y[?];
    dGplasma = -p.k_Gplasma_to_Gbeta * Gplasma
    dGbeta = p.k_Gplasma_to_Gbeta * Gplasma 
    dInsulin = 0

    k_insulin_secretion = p.k_baseline_insulin_secretion
    if Gplasma > NORMAL_BLOOD_SUGAR_HIGH:
        INSULIN_ACTIVE = True
        GLUCAGON_ACTIVE = False
        dinsulin += k_insulin_secretion

    # dydt[?] = dGplasma
    # dydt[?] = dGbeta
    # dydt[?] = dInsulin

    return dydt

def alpha_cell(t, y, p):
    dydt = np.zeros((n, 1))

    # NOTE : insulin activity has 90% chance to inhibit alpha cell activity
    switch = random.random()
    if INSULIN_ACTIVE and switch > 0.1:
        return dydt

    # Gplasma = y[?];
    dGplasma = -p.k_Gplasma_to_Galpha * Gplasma
    dGalpha = p.k_Gplasma_to_Galpha * Gplasma 
    dGlucagon = 0

    k_glucagon_secretion = p.k_baseline_glucagon_secretion
    if Gplasma < NORMAL_BLOOD_SUGAR_LOW:
        GLUCAGON_ACTIVE = True
        INSULIN_ACTIVE = False
        dGlucagon += k_glucagon_secretion

    # dydt[?] = dGplasma
    # dydt[?] = dGbeta
    # dydt[?] = dGlucagon

    return dydt

def delta_cell(t, y, p):
    dydt = np.zeros((n, 1))
    return dydt

def gamma_cell(t, y, p):
    dydt = np.zeros((n, 1))
    return dydt
