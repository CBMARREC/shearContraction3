__author__ = 'Georgios Grekas (grekas.g@gmail.com)'
results_path = 'results/'

from init_parameters import *
import os
import math


# --------- IMPORTANT -------------------------------------

# Body starts with an initial_shear, will be tilted so that the top advances in small steps of length steps_length until 
# it reaches max_shear. From there, it wil "deshear" until reaching final_shear in case hysteresis is intended


# -------- Shear values -----------------------------

initial_shear = 0.04 # this is tan theta (theata is the angle between the lateral side at the beginning and the same side at the end)
max_shear = 0.25 # TOP side will be displaced a (max_shearx100)% from the initial position
final_shear = 0.04 # Final shear after deshearing the body

step_length = 0.0001 # Shear right and deshear left will happen in small steps. Shear for step N+1 has shear of step N as initial guess

print_Niter = 1

resolution = 30

uc0 = -0.0 # uniform radial displacement, - gives contraction, + gives expansion 
###u0 can also be a vector, i.e. [a, b] for some a and b,  or an expression.

my_k = 30 # Value of the constant k

my_model = 'linear_spring' # 'PNIPAAm_particles'

# choose if the centers are free to move or not. Cells can move only if one has called the methodproblem.set_k
solver_t = 'free_centers' # 'fixed_centers' 

el_order = 1

# ---------- Management of number of graphics -----------

number_images = 7
image_frequency = math.ceil((max_shear - final_shear)/(step_length*number_images))

# ---------- end of graphic management ---------

# --------- Create domain and mesh START ---------------------

# Square side

square_side = 10

# Define the body boundaries so later we can apply boaundary conditions over them.
bottom =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0)
top = CompiledSubDomain("near(x[1], side) && on_boundary", side = square_side)#

left = CompiledSubDomain("near(x[0], side) && on_boundary", side = 0)#
right = CompiledSubDomain("near(x[0], side) && on_boundary", side = square_side)


# The cell

#c_x, c_y, rho = 0.5, 0.5, 0.3 # center and radius of a big cirlce (The ECM)  NO EN MI CASO
c1_x, c1_y, i_rho1 = 5., 5., 0.5    # center and radius of a cell

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


print('---------------------------------')
print('------------- PASO 1 ------------')
print('---------------------------------')


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



# ------------- Open necessary files ---------------------------

energies = np.load("results/vtkFiles/saved_functions/energy_i.npy", "r") #opens the file where energies are stored

energyVsShear = open("results/vtkFiles/saved_functions/EnergyVsShear.csv", "w") 

last_line1st = energies[-2] # last energy
energyVsShear.write(str(initial_shear)+','+str(last_line1st)+'\n')

# -------------Ends files management -----------------

# ---------------- solving Ends -----------------------------------------

count = 1

# ------------------- loop begins (shear) ---------------------



while (initial_shear < max_shear):

    u_0 = u
    count = count + 1
    shear = initial_shear + step_length
    initial_shear = shear


    print('---------------------------------')
    print '------------- PASO', count, '------------'
    print '        Shear', shear
    print('---------------------------------')

    # Define Dirichlet boundary (y = 0 or y = 1)
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

    domain = RectangularDomain(0, 0, square_side, square_side, [ [b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False)

    domain.remove_subdomain(cell_domain)

    domain.create_mesh(resolution) # create a mesh over the domain for the given resolution


    problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                        12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
                                         domain, el_order = el_order)

    problem.set_k(my_k, cell_model = my_model) # accounts for cell response  (you can exclude this line)


    solver = m_ncgSolver(problem, res_path = results_path, solver = solver_t, print_Niter = print_Niter) # let save_run = False

    # initialize the displacement vector using a the previous result
    solver.init_from_function(u_0)

    u = solver.solve()

    try:
	    solver.save_all_function_info()
    except:
	    print('hdf5 is not supported')

    if count % image_frequency == 0:
        solver.plot_results()
        file_number=str(shear)
        os.rename(r'results/vtkFiles/u.pvd',r"results/vtkFiles/u_c_"+file_number+".pvd")
        os.rename(r'results/vtkFiles/detF.pvd',r"results/vtkFiles/detF_c_"+file_number+".pvd")
        os.rename(r'results/vtkFiles/dens.pvd',r"results/vtkFiles/dens_c_"+file_number+".pvd")
        os.rename(r'results/vtkFiles/u000000.vtu',r"results/vtkFiles/u_c_"+file_number+".vtu")
        os.rename(r'results/vtkFiles/detF000000.vtu',r"results/vtkFiles/detF_c_"+file_number+".vtu")
        os.rename(r'results/vtkFiles/dens000000.vtu',r"results/vtkFiles/dens_c_"+file_number+".vtu")

    energies = np.load("results/vtkFiles/saved_functions/energy_i.npy", "r") #opens the file where energies are stored

    last_line1st = energies[-2] # last energy
    energyVsShear.write(str(initial_shear)+','+str(last_line1st)+'\n')


# ------------------- loop ends -----------------------


count = 0

# ------------------- 2nd loop begins (deshear) ---------------------

while (final_shear < initial_shear):

    u_0 = u
    count = count + 1
    shear = initial_shear - step_length
    initial_shear = shear


    print('---------------------------------')
    print '------------- PASO', count, 'B  ------------'
    print '        Shear', shear
    print('---------------------------------')

    # Define Dirichlet boundary (y = 0 or y = 1)
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

    domain = RectangularDomain(0, 0, square_side, square_side, [ [b,bottom], [l,left], [r,right], [t, top]], isMeshUniform = False)

    domain.remove_subdomain(cell_domain )

    domain.create_mesh(resolution) # create a mesh over the domain for the given resolution


    problem = mEnergyMinimizationProblem('1.0/96.0 * (5*Ic**3 - 9*Ic**2 - \
                                        12*Ic*J**2 + 12*J**2 + 8) + exp( 80*(0.22 - J) )',
                                         domain, el_order = el_order)

    problem.set_k(my_k, cell_model = my_model) # accounts for cell response  (you can exclude this line)


    solver = m_ncgSolver(problem, res_path = results_path, solver = solver_t, print_Niter = print_Niter) # let save_run = False

    # initialize the displacement vector using a the previous result
    solver.init_from_function(u_0)

    u = solver.solve()

    try:
	    solver.save_all_function_info()
    except:
	    print('hdf5 is not supported')

    if count % image_frequency == 0:
        solver.plot_results()
        file_number=str(shear)
        os.rename(r'results/vtkFiles/u.pvd',r"results/vtkFiles/u_d_"+file_number+".pvd")
        os.rename(r'results/vtkFiles/detF.pvd',r"results/vtkFiles/detF_d_"+file_number+".pvd")
        os.rename(r'results/vtkFiles/dens.pvd',r"results/vtkFiles/dens_d_"+file_number+".pvd")
        os.rename(r'results/vtkFiles/u000000.vtu',r"results/vtkFiles/u_d_"+file_number+".vtu")
        os.rename(r'results/vtkFiles/detF000000.vtu',r"results/vtkFiles/detF_d_"+file_number+".vtu")
        os.rename(r'results/vtkFiles/dens000000.vtu',r"results/vtkFiles/dens_d_"+file_number+".vtu")

    energies = np.load("results/vtkFiles/saved_functions/energy_i.npy", "r") #opens the file where energies are stored

    last_line1st = energies[-2] # last energy
    energyVsShear.write(str(initial_shear)+','+str(last_line1st)+'\n')


# ------------------- loop ends -----------------------

solver.plot_results()


"""
energyVsShear.close()

solver.plot_results()
try:
	solver.save_all_function_info()
except:
	print('hdf5 is not supported')

"""