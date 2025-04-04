# distutils: language = c
# cython: boundscheck=False, wraparound=False, nonecheck=False

from libc.math cimport exp

cdef double ligand_bound(double L, double R, double kon, double koff):
    return (kon * L * R) / koff

cdef double fractional_occupancy_Kd(double L, double Kd):
    return L / (L + Kd)

cdef double fractional_occupancy_LR(double R, double LR):
    return 1.0 / (1.0 + (R / LR))

cdef double fractional_inhibitory_effect(double f_occupancy, double k):
    return 1.0 / (1.0 + exp(k * (f_occupancy - 0.5)))

cdef double inhibitory_secretion(double L, double R, double LR, double kon, double koff, double kmax_secretion, double k):
    cdef double Kd = (L * R) / LR
    cdef double fo = L / (L + Kd)
    cdef double f = 1.0 / (1.0 + exp(k * (fo - 0.5)))
    return kmax_secretion * f

cdef tuple ligandODE(double L, double R, double LR, double kon, double koff):
    cdef double dL = -kon * L * R + koff * LR
    cdef double dR = -kon * L * R + koff * LR
    cdef double dLR = kon * L * R - koff * LR
    return dL, dR, dLR
