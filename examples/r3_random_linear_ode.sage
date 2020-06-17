# Import the code from the src folder without installation
# by temporarily adding the src folder to path.
import sys
import pathlib
path = pathlib.Path().absolute().parent
sys.path.append(str(path))

from ode_abnormals import *

# Define variables for the ODE in a polynomial ring.
# The names of the variables are not important as they will be replaced.
PR = PolynomialRing(QQ, ['x_%d' % (i + 1) for i in range(3)])
ode_variables = PR.gens()

# Form an ODE with random integer coefficients.
ode = [sum(ZZ.random_element() * xk for xk in ode_variables)
       for xk in ode_variables]

# Enable verbose output to get more information about the computation.
set_verbose(1)

# Compute all admissible covector for trajectories of the ODE reaching zero.
# solution_type must be one of "all", "any", or "nice"
solution_type = "nice"
covec, auxdata = ode_search(ode, ode_variables, solution_type)

# Print the ODE.
print("\nODE:")
for xk, dk in zip(ode_variables, ode):
    print("d%s/dt = %s" % (xk, dk))

print("\nA priori: trajectory is abnormal in step %d." % auxdata.step_max)
print("A posteriori: abnormal in step %d." % auxdata.step)

# Print the covector in a text data format.
# The 'm' parameter changes the way the coefficients are labeled.
# A Hall word such as 2221112 would factor as
#   m=3: (2)^3(1)(112)
#   m=2: (2)^3(1)^2(12)
#   m=1: (2)^3(1)^3(2)
#   m=0: 2221112
m = auxdata.layer
print("\nAbnormal covector for the linear ODE:\n")
print(covector_txt(covec, m))

PR = auxdata.polynomial_ring
P = abnormal_polynomials_with_covec(PR, covec, auxdata.layer)
print("\n\nAbnormal polynomials of layer %d for the covector:\n" % m)
for X in P:
    if P[X]:
        pf = factor(P[X])
    else:
        pf = P[X]
    print("P_{%s} = %s" % (indexword(X.leading_support()), pf))
