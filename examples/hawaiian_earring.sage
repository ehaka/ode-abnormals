# Import the code from the src folder without installation
# by temporarily adding the src folder to path.
import sys
import pathlib
path = pathlib.Path().absolute().parent
sys.path.append(str(path))

from ode_abnormals import *

# Define variables for the ODE in a polynomial ring.
# The names of the variables are not important as they will be replaced.
PR = PolynomialRing(QQ, ['x', 'y'])
ode_variables = PR.gens()
x, y = ode_variables

# Form the ODE. The equations are given in the order
# of the variables of the polynomial ring.
ode = [
    x * x - y * y,  # dx/dt = x^2-y^2
    2 * x * y  # dy/dt = 2*x*y
]

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
print("\nA priori: trajectory is abnormal in step %d by coefficients of the series" % auxdata.step_max)
print("%s = %s" % (auxdata.generating_function, auxdata.series))
print("A posteriori: abnormal in step %d." % auxdata.step)

# Print the simplified covector in latex format.
# The 'm' parameter changes the way the coefficients are labeled.
# A Hall word such as 2221112 would factor as
#   m=3: (2)^3(1)(112)
#   m=2: (2)^3(1)^2(12)
#   m=1: (2)^3(1)^3(2)
#   m=0: 2221112
m = auxdata.layer
print("\nAbnormal covector for the Hawaiian earring:\n")
print(covector_latex(covec, m))

# Print the abnormal polynomials of the relevant layer (i.e. layer 4).
P = abnormal_polynomials_with_covec(auxdata.polynomial_ring, covec, m)
print("\n\nAbnormal polynomials of layer %d:\n" % m)
for X in P:
    print("P_{%s} = %s" % (indexword(X.leading_support()), factor(P[X])))
