from enum import IntEnum, auto

class Index(IntEnum):
    """
    Index enumeration for the system's state vector.

    Each attribute represents the index of a specific biological quantity 
    (plasma concentrations, intracellular metabolites, gut contents, etc.)
    within the full state vector used in ODE system simulations.

    Naming convention:
        - plasma_*: variables in blood plasma
        - subq_*: variables in subcutaneous adipose tissue
        - vsc_*: variables in visceral adipose tissue
        - muscle_*: variables in skeletal muscle
        - gut_*: nutrients in the gut
        - micellar_*, membrane_*, cytosol_*: stages of fatty acid absorption and storage in gut

    Auto-incremented indices ensure that each variable is assigned a unique position 
    in the vector.

    Used primarily for mapping between biological compartments and their corresponding 
    positions in the simulation state array.
    """
    plasma_glucose=0
    plasma_fructose=auto()
    plasma_fattyacid=auto()
    plasma_aminoacid=auto()
    plasma_lactate=auto()
    plasma_insulin=auto()
    plasma_glucagon=auto()
    plasma_somatostatin=auto()

    subq_glucose=auto()
    subq_fattyacid=auto()
    subq_aminoacid=auto()
    subq_insulin=auto()
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

    liver_glucose=auto()
    liver_fructose=auto()
    liver_fattyacid=auto()
    liver_aminoacid=auto()
    liver_TAG = auto()
    liver_pyruvate=auto()
    liver_G6P=auto()
    liver_ROS=auto()
    liver_ACoA=auto()

def get_plasma_indices():
    return (
        Index.plasma_glucose,
        Index.plasma_insulin,
        Index.plasma_fattyacid,
        Index.plasma_aminoacid,
        Index.plasma_lactate,
    )

def get_plasma_names():
    return (
        Index.plasma_glucose.name,
        Index.plasma_insulin.name,
        Index.plasma_fattyacid.name,
        Index.plasma_aminoacid.name,
        Index.plasma_lactate.name,
    )

def get_subq_indices():
    return (
        Index.subq_glucose,
        Index.subq_insulin,
        Index.subq_fattyacid,
        Index.subq_aminoacid,
        Index.subq_G6P,
        Index.subq_TAG,
        Index.subq_ACoA,
        Index.subq_pyruvate,
        Index.subq_ROS,
    )

def get_subq_names():
    return (
        Index.subq_glucose.name,
        Index.subq_insulin.name,
        Index.subq_fattyacid.name,
        Index.subq_aminoacid.name,
        Index.subq_G6P.name,
        Index.subq_TAG.name,
        Index.subq_ACoA.name,
        Index.subq_pyruvate.name,
        Index.subq_ROS.name,
    )

def get_vsc_indices():
    return (
        Index.vsc_glucose,
        Index.vsc_insulin,
        Index.vsc_fattyacid,
        Index.vsc_aminoacid,
        Index.vsc_G6P,
        Index.vsc_TAG,
        Index.vsc_ACoA,
        Index.vsc_pyruvate,
        Index.vsc_ROS,
    )
def get_vsc_names():
    return (
        Index.vsc_glucose.name,
        Index.vsc_insulin.name,
        Index.vsc_fattyacid.name,
        Index.vsc_aminoacid.name,
        Index.vsc_G6P.name,
        Index.vsc_TAG.name,
        Index.vsc_ACoA.name,
        Index.vsc_pyruvate.name,
        Index.vsc_ROS.name,
    )

def get_muscle_indices():
    return (
        Index.muscle_glucose,
        Index.muscle_insulin,
        Index.muscle_fattyacid,
        Index.muscle_aminoacid,
        Index.muscle_G6P,
        Index.muscle_glycogen,
        Index.muscle_pyruvate,
        Index.muscle_ACoA,
        Index.muscle_NAD,
        Index.muscle_NADH,
        Index.muscle_FAD,
        Index.muscle_FADH2,
        Index.muscle_ROS,
        Index.muscle_ATP,
        Index.muscle_lactate,
    )

def get_muscle_names():
    return (
        Index.muscle_glucose.name,
        Index.muscle_insulin.name,
        Index.muscle_fattyacid.name,
        Index.muscle_aminoacid.name,
        Index.muscle_G6P.name,
        Index.muscle_glycogen.name,
        Index.muscle_pyruvate.name,
        Index.muscle_ACoA.name,
        Index.muscle_NAD.name,
        Index.muscle_NADH.name,
        Index.muscle_FAD.name,
        Index.muscle_FADH2.name,
        Index.muscle_ROS.name,
        Index.muscle_ATP.name,
        Index.muscle_lactate.name,
    )

def get_gut_indices():
    return (
        Index.gut_glucose,
        Index.gut_fructose,
        Index.micellar_fattyacid,
        Index.membrane_fattyacid,
        Index.cytosol_fattyacid,
        Index.cytosol_TAG,
    )
def get_gut_names():
    return (
        Index.gut_glucose.name,
        Index.gut_fructose.name,
        Index.micellar_fattyacid.name,
        Index.membrane_fattyacid.name,
        Index.cytosol_fattyacid.name,
        Index.cytosol_TAG.name,
    )
