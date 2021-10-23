__author__ = "Georgios Grekas (grekas.g@gmail.com)"
results_path = "results/"

#############################
#### Initial Parameters #####
#############################

from dolfin import *
from mClasses.DomainCs import *
from mClasses.ProblemCs import *
from mClasses.SolverCs import *

# Optimization options for the form compiler
parameters["mesh_partitioner"] = "SCOTCH"  # "ParMETIS" #
parameters["partitioning_approach"] = "PARTITION"
#
# # parameters["form_compiler"]["precision"] = 1000
parameters["allow_extrapolation"] = True
# Make mesh ghosted for evaluation of DG terms
parameters["ghost_mode"] = "shared_facet"  # Other options: 'shared_vertex', 'none'
form_compiler_parameters = {"quadrature_degree": 2}

parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"
ffc_options = {
    "optimize": True,
    "eliminate_zeros": True,
    "precompute_basis_const": True,
    "precompute_ip_const": True,
}
# --------- IMPORTANT -------------------------------------

# Body starts with an initial_shear, will be tilted so that the top advances in small steps of length steps_length until
# it reaches max_shear. From there, it wil decompress until reaching final_shear in case hysteresis is intended

initial_shear = 0.05
max_shear = 0.2  # TOP side will be displaced a 100% from the initial position

step_lenght = 0.005

print_Niter = 1

mesh_Resolution = 6

# --------- Create domain and mesh START ---------------------

# Le cuadrado
bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side=0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side=1)

left = CompiledSubDomain("near(x[0], side) && on_boundary", side=0)
right = CompiledSubDomain("near(x[0], side) && on_boundary", side=1)

inside = CompiledSubDomain("near(x[0], side) && on_boundary", side=1)

# The cell

# c_x, c_y, rho = 0.5, 0.5, 0.3 # center and radius of a big cirlce (The ECM)  NO EN MI CASO
c1_x, c1_y, i_rho1 = 0.5, 0.5, 0.3  # center and radius of a cell
# c2_x, c2_y, i_rho2 = -4, 0., 2.   # center and radius of the second cell


cell_domain = CircularDomain(
    i_rho1, c1_x, c1_y, 0
)  # creates a circle with radius rho and center at (c_x, c_y), 0 indicates that
# the outer boundary is fixed
# domain =CircularDomain(rho, c_x, c_y) # here the outer boundary is free


# Define Dirichlet boundary (x = 0 or x = 1)
b = Expression(("0.0", "0.0"), pr=initial_shear, degree=2)
t = Expression(("1*pr", "0.0"), pr=initial_shear, degree=2)

l = Expression(("pr*x[1]", "0.0"), pr=initial_shear, degree=2)

r = Expression(("pr*x[1]", "0.0"), pr=initial_shear, degree=2)


# u0 = -0.5 # uniform radial displacement, - gives contraction, + gives expansion
###u0 can also be a vector, i.e. [a, b] for some a and b,  or an expression.

domain = RectangularDomain(
    0, 0, 1, 1, [[b, bottom], [l, left], [r, right], [t, top]], isMeshUniform=False
)

domain.remove_subdomain(CircularDomain(i_rho1, c1_x, c1_y, 0))
# isMeshUnifor default value is "False". When "True" the area of all triangles in the original mesh are equal. Otherwise triangles ares are as better fit.
# !!! When there are irregular shapes inside (ex: circles for cells) this value has to be False so triangles can addapt to the irregular form.

resolution = (
    mesh_Resolution  # determine the mesh resolution, higher resolution gives finer mesh
)
domain.create_mesh(resolution)  # create a mesh over the domain for the given resolution
mesh = domain.get_mesh()
# --------- Create domain and mesh  END ---------------------

# -----------------Define problem type START --------------------------------
# ---------- i.e. mEnergyMinimization or mNonlinearVariational --------------
el_order = 1  # the elements order
problem = mEnergyMinimizationProblem(
    "1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                    12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )",
    domain,
    el_order=el_order,
)
# problem.set_k(1, cell_model ='linear_spring') # accounts for cell response
# problem.set_k(1, cell_model ='PNIPAAm_particles') # accounts for cell response
# -----------------Define problem type END --------------------------------

# ---------------- Define solver type--------------------------------------

# choose if the centers are free to move or not. Cells can move only if one has called the method problem.set_k
solver_t = "fixed_centers"  #'free_centers'
solver = m_ncgSolver(
    problem, res_path=results_path, solver=solver_t, print_Niter=print_Niter
)  # let save_run =False,

print("---------------------------------")
print("------------- PASO 1 ------------")
print("---------------------------------")

solver.initialization(
    "Polyconvex"
)  # initialize the displacement vector using a polyconvex functional
solver.add_disturbances()  # If the program finds a SADdle point or local minimum, disturbances will take the solution out of it and the solution can find the global min.
u = solver.solve()

try:
    solver.save_all_function_info()
except:
    print("hdf5 is not supported")

# ------------- Open necessary files ---------------------------

energies = np.load(
    "results/vtkFiles/saved_functions/energy_i.npy", "r"
)  # opens the file where energies are stored

energyVsShear = open("results/vtkFiles/saved_functions/EnergyVsShear.csv", "w")

# last_line2nd = energies[-3] # second to the last energy
last_line1st = energies[-2]  # last energy
# energyVsCompression.write(str(initial_compression)+','+str(last_line2nd)+'\n')
energyVsShear.write(str(initial_shear) + "," + str(last_line1st) + "\n")

# -------------Ends files management -----------------

# ---------------- solving Ends -----------------------------------------

count = 1

# ------------------- loop begins ---------------------


while initial_shear < max_shear:

    u_0 = u
    count = count + 1
    shear = initial_shear + step_lenght
    initial_shear = shear

    print("---------------------------------")
    print(("------------- PASO", count, "------------"))
    print(("        Shear", shear))
    print("---------------------------------")

    bottom = CompiledSubDomain("near(x[1], side) && on_boundary", side=0)
    top = CompiledSubDomain("near(x[1], side) && on_boundary", side=1)  #

    left = CompiledSubDomain("near(x[0], side) && on_boundary", side=0)  #
    right = CompiledSubDomain("near(x[0], side) && on_boundary", side=1)

    # Define Dirichlet boundary (y = 0 or y = 1)
    b = Expression(("0.0", "0.0"), pr=initial_shear, degree=2)
    t = Expression(("1*pr", "0.0"), pr=initial_shear, degree=2)

    l = Expression(("pr*x[1]", "0.0"), pr=initial_shear, degree=2)

    r = Expression(("pr*x[1]", "0.0"), pr=initial_shear, degree=2)

    # domain = RectangularDomain(0, 0, 1, 1, [ [b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False)

    resolution = mesh_Resolution  # determine the mesh resolution, higher resolution gives finer mesh

    # domain.create_mesh(resolution) # create a mesh over the domain for the given resolution
    domain.create_mesh(
        resolution
    )  # create a mesh over the domain for the given resolution
    mesh = domain.get_mesh()

    problem = mEnergyMinimizationProblem(
        "1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                        12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )",
        domain,
        el_order=el_order,
    )

    # choose if the centers are free to move or not. Cells can move only if one has called the method problem.set_k
    solver_t = "fixed_centers"  #'free_centers'

    solver = m_ncgSolver(
        problem, res_path=results_path, solver=solver_t, print_Niter=print_Niter
    )  # let save_run =False

    # initialize the displacement vector using a the previous result
    solver.init_from_function(u_0)

    u = solver.solve()

    try:
        solver.save_all_function_info()
    except:
        print("hdf5 is not supported")

    energies = np.load("results/vtkFiles/saved_functions/energy_i.npy", "r")
    #   last_line2nd = energies[-3] # second to the last energy
    last_line1st = energies[-2]  # last energy
    #   energyVsCompression.write(str(initial_compression)+','+str(last_line2nd)+'\n')
    energyVsShear.write(str(initial_shear) + "," + str(last_line1st) + "\n")


# ------------------- loop ends -----------------------


energyVsShear.close()


solver.plot_results()
try:
    solver.save_all_function_info()
except:
    print("hdf5 is not supported")
