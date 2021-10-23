# This program takes the EnergyVsCompression and 1) approximates with B-splines and
# 2) calculates the derivative also with B-splines
import matplotlib
import numpy as np
from numpy import ndarray

import matplotlib.pyplot as plt
from scipy import interpolate

# matplotlib.use("Agg")

energies_vs_strains = np.genfromtxt(
    "EnergyVsShear.csv", delimiter=","
)  # type: ndarray # Imports the data in the file EnergyVsCompression.csv

strains = energies_vs_strains[
          0:-1, 0
          ]  # From the previous file takes only the first column which is the strains
energies = energies_vs_strains[
           0:-1, 1
           ]  # From the previous file takes only the second column, which is the energies

max_index = np.argmax(
    strains
)  # Finds the position in the list strains with the higher strain.
# The goal is to separate the list in two parts: the part of shear and the part of unshear.

strains_shear = strains[
                      0:max_index
                ]  # It takes only the strains corresponding to the shear
energies_shear = energies[
                       0:max_index
                 ]  # It takes only the energies corresponding to the shear part

strains_unshear = np.flipud(
    strains[max_index + 1: -1]
)  # It takes only the strains corresponding to the unshear part
energies_unshear = np.flipud(
    energies[max_index + 1: -1]
)  # It takes only the energies corresponding to the unshear part

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


# =============== END OF AKIMA PART ============

#

# ============= STRESS PRODUCED FROM NUMERICAL DERIVATION OF INTERPOLATED ENERGY ==================

count = 0
count2 = 0

data_auxiliar_shear = np.zeros((stress_Akima_shear.size, 2), dtype=float)
data_auxiliar_unshear = np.zeros(
    (stress_Akima_unshear.size, 2), dtype=float
)

for x in stress_Akima_shear:
    data_auxiliar_shear[count, 0] = strains_shear[count]
    data_auxiliar_shear[count, 1] = stress_Akima_shear[count]
    count += 1  # esto es count = count +1

for x in stress_Akima_unshear:
    data_auxiliar_unshear[count2, 0] = strains_unshear[count2]
    data_auxiliar_unshear[count2, 1] = stress_Akima_unshear[count2]
    count2 += 1  # esto es count = count +1

# """""" Files that store the points for stress vs strain """"""
np.savetxt(
    "stress_Akima_shear.csv",
    data_auxiliar_shear,
    delimiter=",",
    fmt=["%.11f", "%.11f"],
)
np.savetxt(
    "stress_Akima_unshear.csv",
    data_auxiliar_unshear,
    delimiter=",",
    fmt=["%.11f", "%.11f"],
)

# """""" Plots """"""
# plt.plot(strains_shear[3500:3700], stress_Akima_shear[3500:3700], "go")
# plt.savefig("stress_Akima_zoom", dpi=900)

# ============ END OF NUM STRESS FROM INTERPL ENERGY =========================================

#

# ========================= STRESS FROM INTERPOLATION =======================================


# """""" The interpolator method derivative gives a interpolated function for the nu-esima derivative """"""""
# stress_interpolator_compression = energies_interpolator_shear.derivative(nu=1)
# stress_interpolate = energies_interpolator_shear.derivative(1)

# plt.plot(strains_shear, stress_interpolate(strains_shear))
# plt.savefig("stress_interpolated")

# plt.plot(strains_shear[2000:5000], stress_interpolate(strains_shear[2000:5000]))
# plt.savefig("stress_interpolated_zoom")

# plt.plot(
#     strains_shear[3500:3700], stress_interpolate(strains_shear[3500:3700])
# )
# plt.savefig("stress_interpolated_comparison2.eps", dpi=900)

# plt.plot(strains_shear, energies_interpolator_shear(strains_shear))
# plt.savefig("stress_interp.pdf")
