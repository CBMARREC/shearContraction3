__author__ = 'Georgios Grekas (grekas.g@gmail.com)'
results_path = 'results2/'

from init_parameters import *
import os
import math


# --------- IMPORTANT -------------------------------------

# Body suffers shear deformation = initial_shear in only one step all at once


# -------- Shear values -----------------------------

initial_shear = 0.25 # this is tan theta (theata is the angle between the lateral side at the beginning and the same side at the 

print_Niter = 1

resolution = 100

square_side = 10

uc0 = -0.50 # uniform radial displacement, - gives contraction, + gives expansion 
###u0 can also be a vector, i.e. [a, b] for some a and b,  or an expression.

my_k = 30 # Value of the constant k

my_model = 'linear_spring' # 'PNIPAAm_particles'

# choose if the centers are free to move or not. Cells can move only if one has called the methodproblem.set_k
solver_t = 'free_centers' # 'fixed_centers' 

el_order = 1



# --------- Create domain and mesh START ---------------------

# Define the body boundaries so later we can apply boaundary conditions over them.
bottom =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = square_side)#

left = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)#
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = square_side)
# inside = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1)

# The cell

#c_x, c_y, rho = 0.5, 0.5, 0.3 # center and radius of a big cirlce (The ECM)  NO EN MI CASO
c1_x, c1_y, i_rho1 = 5., 5., 1.    # center and radius of a cell

# cell_domain = CircularDomain(i_rho1, c1_x, c1_y) # creates a circle with radius rho and center at (c_x, c_y), 0 indicates that
# the outer boundary is fixed
# cell_domain =CircularDomain(i_rho1, c1_x, c1_y) # here the outer boundary is free
cell_domain = CircularDomain(i_rho1, c1_x, c1_y, uc0) # uc0 is the restriction over the cell boundary, for example a contraction


# Define Dirichlet boundary (x = 0 or x = 1)
b = Expression(("0.0",
                "0.0"),pr = initial_shear, degree = 2)
t = Expression(("10*pr",
                "0.0"),
                pr = initial_shear, degree = 2)

l = Expression(("pr*x[1]",
                "0.0"),
                pr = initial_shear, degree = 2)

r = Expression(("pr*x[1]",
                "0.0"),
                pr = initial_shear, degree = 2)


domain = RectangularDomain(0, 0, square_side, square_side, [ [b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False) #isMeshUnifor default value is "False". When "True" the area of all triangles in the original mesh are equal. Otherwise triangles ares are as better fit.
domain.remove_subdomain(cell_domain)
# !!! When there are irregular shapes inside (ex: circles for cells) this value has to be False so triangles can addapt to the irregular form.

domain.create_mesh(resolution) # create a mesh over the domain for the given resolution
mesh = domain.get_mesh()

# --------- Create domain and mesh  END ---------------------



# -----------------Define problem type START --------------------------------
# ---------- i.e. mEnergyMinimization or mNonlinearVariational --------------

problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                    12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
                                     domain, el_order = el_order)
 
problem.set_k(my_k, cell_model = my_model) # accounts for cell response  

solver = m_ncgSolver(problem, res_path = results_path, solver = solver_t, print_Niter = print_Niter) # let save_run = False, this has to be in the while loop every time

# -----------------Define problem type END --------------------------------


solver.initialization('Polyconvex') # initialize the displacement vector using a polyconvex functional
solver.add_disturbances() # If the program finds a SADdle point or local minimum, disturbances will take the solution out of it and the solution can find the global min.
u = solver.solve()

solver.plot_results()

try:
	solver.save_all_function_info()
except:
	print('hdf5 is not supported')


# ---------------- solving Ends -----------------------------------------



solver.plot_results()


"""
energyVsShear.close()

solver.plot_results()
try:
	solver.save_all_function_info()
except:
	print('hdf5 is not supported')

"""