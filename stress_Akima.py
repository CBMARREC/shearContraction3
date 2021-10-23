# This program takes the EnergyVsShear and
# 1) approximates with B-splines and
# 2) calculates the derivative also with B-splines
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate

matplotlib.use("Agg")


# Imports the data in the file EnergyVsShear.csv
energies_vs_strains = np.genfromtxt("EnergyVsShear.csv", delimiter=",")
# From the previous file takes only the first column which is the strains
strains = energies_vs_strains[0:-1, 0]
# From the previous file takes only the second column, which is the energies
energies = energies_vs_strains[0:-1, 1]
# Finds the position in the list strains with the higher strain.
# The goal is to separate the list in two parts: the part of shear and the part of unshear.
max_index = np.argmax(strains)

# It takes only the strains corresponding to the shear
strains_shear = strains[0:max_index]
# It takes only the energies corresponding to the shear part
energies_shear = energies[0:max_index]
# It takes only the strains corresponding to the unshear
strains_unshear = np.flipud(strains[max_index + 1: -1])
# It takes only the energies corresponding to the unshear part
energies_unshear = np.flipud(energies[max_index + 1: -1])


# ================ AKIMA PART =================

# """""" Constructing two objects to interpolate data from imported files """"""""
energies_interpolator_shear = interpolate.Akima1DInterpolator(
    strains_shear, energies_shear, axis=0
)
energies_interpolator_unshear = interpolate.Akima1DInterpolator(
    strains_unshear, energies_unshear, axis=0
)

# """"" The interpolator method __call__ allows to compute the derivative nu-esima of the interpolate object """""""
stress_Akima_shear = energies_interpolator_shear.__call__(
    strains_shear, nu=1, extrapolate=None
)
stress_Akima_unshear = energies_interpolator_unshear.__call__(
    strains_unshear, nu=1, extrapolate=None
)

# """""" The interpolator method derivative gives a interpolated function for the nu-esima derivative """"""""
stress_interpolator_shear = energies_interpolator_shear.derivative(nu=1)


# =============== END OF AKIMA PART ============

#

# ============= STRESS PRODUCED FROM NUMERICAL DERIVATION OF INTERPOLATED ENERGY ==================

def data_auxiliar_matrix(stress_akima, file_name):
    data_auxiliar = np.zeros((stress_akima.size, 2), dtype=float)
    for count in range(len(stress_akima)):
        data_auxiliar[count, 0] = strains[count]
        data_auxiliar[count, 1] = stress_akima[count]

    np.savetxt(
        file_name,
        data_auxiliar,
        delimiter=",",
        fmt=["%.11f", "%.11f"],
    )


# """""" Files that store the points for stress vs strain """"""
data_auxiliar_matrix(stress_Akima_shear, "stress_Akima_shear.csv")
data_auxiliar_matrix(stress_Akima_unshear, "stress_Akima_unshear.csv")

# """""" Plots """"""
plt.plot(strains_shear, stress_Akima_shear, "go")
plt.savefig("stress_Akima_zoom", dpi=900)


# ============ END OF NUM STRESS FROM INTERPL ENERGY =========================================

#

# ========================= STRESS FROM INTERPOLATION =======================================

stress_interpolate = energies_interpolator_shear.derivative(1)

# plt.plot(strains_shear, stress_interpolate(strains_shear))
# plt.savefig("stress_interpolated")

# plt.plot(strains_shear[2000:5000], stress_interpolate(strains_shear[2000:5000]))
# plt.savefig("stress_interpolated_zoom")

plt.plot(
    strains_shear, stress_interpolate(strains_shear)
)
plt.savefig("stress_interpolated_comparison2.eps", dpi=900)

# plt.plot(strains_shear, energies_interpolator_shear(strains_shear))
# plt.savefig("stress_interp.pdf")
