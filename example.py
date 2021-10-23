__author__ = 'Georgios Grekas (grekas.g@gmail.com)'
results_path = 'results/'

from init_parameters import *

el_order = 1 # the elements order 

# --------- Create domain and mesh START ---------------------
c1_x, c1_y, i_rho1 = 0.5, 0.5, 0.24    # center and radius of a cell

uc0 = -0.10 # uniform radial displacement, - gives contraction, + gives expansion 
###u0 can also be a vector, i.e. [a, b] for some a and b,  or dolfin.Expression.
initial_shear = 0.05

# Le cuadrado
bottom =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = 1)#

left = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)#
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = 1)

# Define Dirichlet boundary (x = 0 or x = 1)
b = Expression(("0.0",
                "0.0"),pr = initial_shear, degree = 2)
t = Expression(("1*pr",
                "0.0"),
                pr = initial_shear, degree = 2)

l = Expression(("pr*x[1]",
                "0.0"),
                pr = initial_shear, degree = 2)

r = Expression(("pr*x[1]",
                "0.0"),
                pr = initial_shear, degree = 2)


domain = RectangularDomain(0, 0, 1, 1, [[b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False)
# domain =CircularDomain(rho, c_x, c_y, 0) # creates a circle with radius rho and center at (c_x, c_y), 0 indicates that
# the outer boundary is fixed
# domain =CircularDomain(rho, c_x, c_y) # here the outer boundary is free

domain.remove_subdomain( CircularDomain(i_rho1, c1_x, c1_y, uc0) ) # from the previous defined domain remove a cirlce 
#                                                                   of radius i_rho1 with center (c1_x, c1_y),
#                                                                   u0 is the boundary condition on the cell

resolution = 6 # determine the mesh resolution, higher resolution gives finer mesh
domain.create_mesh(resolution) # create a mesh over the domain for the given resolution

# --------- Create domain and mesh  END ---------------------


# -----------------Define problem type START --------------------------------
# ---------- i.e. mEnergyMinimization or mNonlinearVariational --------------
problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                    12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
                                     domain, el_order=el_order)
# problem.set_k(1, cell_model ='linear_spring') # accounts for cell response  (you can exclude this line)
# or 
problem.set_k(10, cell_model ='PNIPAAm_particles') # accounts for cell response  (you can exclude this line)
# -----------------Define problem type END --------------------------------




# ---------------- Define solver type--------------------------------------

# choose if the centers are free to move or not. Cells can move only if one has called the method problem.set_k
solver_t = 'free_centers'# 'fixed_centers'#

solver = m_ncgSolver(problem, res_path = results_path, solver=solver_t)
solver.initialization('Polyconvex') # initialize the displacement vector using a polyconvex functional

u = solver.solve()
# ---------------- solving Ends -----------------------------------------




solver.plot_results()
try:
	solver.save_all_function_info()
except:
	print('hdf5 is not supported')


