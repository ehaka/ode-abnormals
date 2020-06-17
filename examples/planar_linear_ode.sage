from sage.all import *

# Import the code from the src folder without installation
# by temporarily adding the src folder to path.
import sys
import pathlib
path = pathlib.Path().absolute().parent
sys.path.append(str(path))

from ode_abnormals import *

# Define free parameters in a fraction field.
CR = PolynomialRing(ZZ, ['a', 'b', 'c', 'd'])
FR = CR.fraction_field()
a11, a12, a21, a22 = FR.gens()

# Define variables for the ODE in a polynomial ring over the fraction field.
# The names of the variables are not important as they will be replaced.
PR = PolynomialRing(FR, ['x', 'y'])
ode_variables = PR.gens()
x, y = ode_variables

# Form the ODE. The equations are given in the order
# of the variables of the polynomial ring.
ode = [
    a11 * x + a12 * y,  # dx/dt
    a21 * x + a22 * y  # dy/dt
]

# Enable verbose output to get more information about the computation.
set_verbose(1)

# Compute a covector for trajectories of the ODE reaching zero.
# solution_type must be one of "all", "any", or "nice"
solution_type = "any"
covec, auxdata = ode_search(ode, ode_variables, solution_type)

# Print the ODE.
print("\nODE:")
for xk, dk in zip(ode_variables, ode):
    print("d%s/dt = %s" % (xk, dk))

# Print the naive abnormality step bound and generating series.
print("\nA priori: trajectory is abnormal in step %d by coefficients of the series" % auxdata.step_max)
print("%s = %s" % (auxdata.generating_function, auxdata.series))
print("A posteriori: abnormal in step %d." % auxdata.step)

# Print the simplified covector as a txt dump.
# The 'm' parameter changes the way the coefficients are labeled.
# A Hall word such as 2221112 would factor as
#   m=3: (2)^3(1)(112)
#   m=2: (2)^3(1)^2(12)
#   m=1: (2)^3(1)^3(2)
#   m=0: 2221112
m = auxdata.layer
print("\nAbnormal covector for a generic homogeneous planar linear ODE:\n")
# Convert the fraction field coefficients to polynomials to enable factoring.
CR = CR.change_ring(QQ)
covec = {X: CR(c.numerator()) / QQ(c.denominator()) for X, c in covec.items()}
covec_txt = covector_txt(covec, m, factor_coeff=True)
print(covec_txt)

# Print the abnormal polynomials of layer m (m=3).
P = abnormal_polynomials_with_covec(auxdata.polynomial_ring, covec, m)
print("\n\nAbnormal polynomials of layer %d:" % m)
for X in P:
    pstr = "\nP_{%s} = " % indexword(X.leading_support())
    padding = len(pstr) - 2
    print(pstr, end="")
    nextpadding = 0
    for mon in P[X].monomials():
        coeff = "(%s)" % (P[X].monomial_coefficient(mon))
        if nextpadding:
            coeff = "+" + coeff
        print("%s%s %s" % (" "*nextpadding, coeff, mon))
        nextpadding = padding

