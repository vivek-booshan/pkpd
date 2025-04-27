from enum import IntEnum, auto

class Index(IntEnum):
    plasma_glucose=0
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
