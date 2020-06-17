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
    -x - y,  # dx/dt = -x -y
    x - y  # dy/dt =  x -y
]

# Enable verbose output to get more information about the computation.
set_verbose(1)

# Compute all admissible covector for trajectories of the ODE reaching zero.
# solution_type must be one of "all", "any", or "nice"
solution_type = "all"
covecs, auxdata = ode_search(ode, ode_variables, solution_type)

# Print the ODE.
print("\nODE:")
for xk, dk in zip(ode_variables, ode):
    print("d%s/dt = %s" % (xk, dk))

print("\nA priori: trajectory is abnormal in step %d by coefficients of the series" % auxdata.step_max)
print("%s = %s" % (auxdata.generating_function, auxdata.series))
print("A posteriori: abnormal in step %d." % auxdata.step)

# Print the covector in a text data format.
# The 'm' parameter changes the way the coefficients are labeled.
# A Hall word such as 2221112 would factor as
#   m=3: (2)^3(1)(112)
#   m=2: (2)^3(1)^2(12)
#   m=1: (2)^3(1)^3(2)
#   m=0: 2221112
m = auxdata.layer
for k, covec in enumerate(covecs):
    print("\nAbnormal covector #%d for the logarithmic spiral:\n" % (k + 1))
    print(covector_txt(covec, m))
print("")

# Alternatively, compute a simplified covector. This attempts to
# decrease the number of nonzero coefficients and their magnitudes.
set_verbose(0)
solution_type = "nice"
covec, auxdata = ode_search(ode, ode_variables, solution_type)

# Print the simplified covector in latex format.
print("\nSimple abnormal covector:\n")
print(covector_latex(covec, 0))

# Print the abnormal polynomials of the relevant layer (i.e. layer 3).
PR = auxdata.polynomial_ring
P = abnormal_polynomials_with_covec(PR, covec, auxdata.layer)
print("\n\nAbnormal polynomials of layer %d for the simple covector:\n" % m)
for X in P:
    print("P_{%s} = %s" % (indexword(X.leading_support()), factor(P[X])))
