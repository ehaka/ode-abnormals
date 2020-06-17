# Warning: the computation of this system goes to rank 3 step 13
# and involves forming and solving a 81360x44383 linear system.
# A sample runtime was 55 minutes 16 seconds
w = 80
print()
print("="*w)
print("="*w)
txt = "    WARNING: COMPUTATION MAY TAKE ROUGHLY AN HOUR    "
l = len(txt)
a = (w - l) // 2
b = w - l - a
print("="*a + txt + "="*b)
print("="*w)
print("="*w)
print()

# Import the code from the src folder without installation
# by temporarily adding the src folder to path.
import sys
import pathlib
path = pathlib.Path().absolute().parent
sys.path.append(str(path))

from ode_abnormals import *

# Define variables for the ODE in a polynomial ring.
# The names of the variables are not important as they will be replaced.
PR = PolynomialRing(QQ, ['x', 'y', 'z'])
ode_variables = PR.gens()
x, y, z = ode_variables

# Form the ODE. The equations are given in the order
# of the variables of the polynomial ring.
lorenz_ode = [
    10 * (y - x),  # dx/dt
    28 * x - y - x * z,  # dy/dt
    x * y - 8 / 3 * z  # dz/dt
]

# Choose an initial point for the trajectory.
initial_point = [1, 0, 0]
# Translate the ODE so that the initial point is the origin.
translation = {xk: xk - ak for xk, ak in zip(ode_variables, initial_point)}
ode = [dx.subs(translation) for dx in lorenz_ode]

# Enable verbose output to get more information about the computation.
set_verbose(1)

# Compute a covector for trajectories of the ODE reaching zero.
# solution_type must be one of "all", "any", or "nice"
# "nice" attempts to decrease the number of nonzero coefficients
# and their magnitudes.
solution_type = "nice"
covec, auxdata = ode_search(ode, ode_variables, solution_type)

# Print the ODE.
print("\nODE:")
for xk, dk in zip(ode_variables, ode):
    print("d%s/dt = %s" % (xk, dk))

# Print the naive abnormality step bound and generating series.
print("\nA priori: trajectory is abnormal in step %d by coefficients of the generating function" % auxdata.step_max)
print(auxdata.generating_function)
print("A posteriori: abnormal in step %d." % auxdata.step)

# Print the simplified covector in latex format.
# The 'm' parameter changes the way the coefficients are labeled.
# A Hall word such as 2221112 would factor as
#   m=3: (2)^3(1)(112)
#   m=2: (2)^3(1)^2(12)
#   m=1: (2)^3(1)^3(2)
#   m=0: 2221112
m = auxdata.layer
print("\nAbnormal covector for the Lorenz butterfly:\n")
print(covector_txt(covec, m))

# Print the abnormal polynomials of the relevant layer (i.e. layer 4).
P = abnormal_polynomials_with_covec(auxdata.polynomial_ring, covec, m)
print("\n\nAbnormal polynomials of layer %d:\n" % m)
for X in P:
    if P[X]:
        Pstr = str(factor(P[X]))
    else:
        Pstr = str(P[X])

    print("P_{%s} = %s" % (indexword(X.leading_support()), Pstr))
