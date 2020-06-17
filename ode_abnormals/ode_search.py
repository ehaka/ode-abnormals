from .abnormals import abnormal_factor_system
from .hall_quotient import NotAFreeLieAlgebra
from .pde import first_integral_pde, integrate_pde

from collections import deque
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from sage.arith.misc import gcd
from sage.calculus.var import var
from sage.functions.log import log
from sage.functions.other import ceil
from sage.misc.defaults import series_precision
from sage.misc.functional import symbolic_prod
from sage.misc.misc import verbose
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.ring import SR
from time import time


class AuxData():
    r"""
    A helper class to contain various data related to computation of abnormal
    covectors.
    """

    def __init__(self, L, m, s, Q, s_max, gf):
        self.lie_algebra = L
        self.rank = len(L.gens())
        self.layer = m
        self.step = s
        self.first_integral = Q
        self.polynomial_ring = Q.parent()
        self.step_max = s_max
        self.generating_function = gf
        self.series = gf.series(gf.variables()[0], s_max + 1)


def _step_upper_bound_low_mem(r, m, q, generating_function):
    r"""
    Low memory implementation of :func:`_step_upper_bound_internal`.

    Significantly slower, but the memory footprint does not significantly
    increase even if the series coefficients need to be computed to very high
    degree terms.
    """
    L = LieAlgebra(ZZ, ['X_%d' % k for k in range(r)]).Lyndon()
    dim_fm = L.graded_dimension(m)

    PR = PolynomialRing(ZZ, 't')
    t = PR.gen()
    a = (1 - dim_fm * (1 - t ** q)) * t ** m
    b = PR.one()
    for k in range(1, m):
        b *= (1 - t ** k) ** L.graded_dimension(k)

    # extract initial coefficients from a symbolic series expansion
    bd = b.degree()
    id = max(a.degree() + 1, bd)
    offset = id - bd
    quot = SR(a / b)
    sym_t = SR(t)
    qs = quot.series(sym_t, id)

    # check if partial sum is positive already within series expansion
    # store the last offset...id terms to start the linear recurrence
    coeffs = deque()
    cumul = ZZ.zero()
    for s in range(id):
        c = ZZ(qs.coefficient(sym_t, s))
        cumul += c
        if s >= offset:
            coeffs.append(c)
        if cumul > 0:
            if generating_function:
                return s, quot
            return s

    # the rest of the coefficients are defined by a recurrence relation
    multipliers = [-b.monomial_coefficient(t ** (bd - k)) for k in range(bd)]
    while cumul <= 0:
        c_next = sum(c * m for c, m in zip(coeffs, multipliers))
        cumul += c_next
        s += 1
        coeffs.append(c_next)
        coeffs.popleft()

    if generating_function:
        return s, quot
    return s


def _step_upper_bound_internal(r, m, q, generating_function):
    r"""
    Common implementation for :func:`step_upper_bound` and :func:`step_upper_bound_specific`
    """
    L = LieAlgebra(QQ, ['X_%d' % k for k in range(r)]).Lyndon()
    dim_fm = L.graded_dimension(m)

    t = var('t')
    pt = symbolic_prod((1 / ((1 - t ** k) ** L.graded_dimension(k))
                        for k in range(1, m)),
                        (1 - dim_fm * (1 - t ** q)) * t ** m)

    increment = series_precision()
    i = 0
    coeffsum = ZZ.zero()
    while True:
        # sum the coefficients of the series for t^i...t^(i+increment-1)
        p_series = pt.series(t, i + increment)
        for s in range(i, i + increment):
            coeffsum += ZZ(p_series.coefficient(t, s))
            if coeffsum > 0:
                verbose("upper bound for abnormality: step %d" % s)
                if generating_function:
                    return s, pt
                return s
        # exponential growth to avoid recomputing the series too much
        i += increment
        increment += increment


def step_upper_bound(r, d, generating_function=False, low_memory=True):
    r"""
    Return an upper bound for the smallest step `s` where trajectories of
    every polynomial ODE system in R^r with polynomials of degree at most `d`
    lift to abnormal curves in the free Carnot group of rank `r` and step `s`.

    INPUT:

    - ``r`` -- the dimension of the ODE system
    - ``d`` -- the maximum degree of polynomials in the system
    - ``generating_function`` -- (default:``False``) a boolean; if ``True``,
      return also the generating function for the series coefficients that
      bound the equation-variable difference
    - ``low_memory`` -- (default:``True``) a boolean; if ``True``,
      uses a low memory footprint method to compute series coefficients.
      Setting to ``False`` and having a well set series_precision can improve
      speed significantly at the cost of memory.
    """
    m = d + 2
    q = d + 1

    if low_memory:
        return _step_upper_bound_low_mem(r, m, q, generating_function)
    return _step_upper_bound_internal(r, m, q, generating_function)


def step_upper_bound_specific(Q, generating_function=False, low_memory=True):
    r"""
    Return a more precise upper bound for the smallest nilpotency step where
    there exists a nontrivial abnormal covector admitting the factor Q.

    INPUT:

    - ``Q`` -- a polynomial to search for as an abnormal factor
    - ``generating_function`` -- (default:``False``) a boolean; if ``True``,
      return also the generating function for the series coefficients that
      bound the equation-variable difference
    - ``low_memory`` -- (default:``True``) a boolean; if ``True``,
      uses a low memory footprint method to compute series coefficients.
      Setting to ``False`` and having a well set series_precision can improve
      speed significantly at the cost of memory.
    """
    PR = Q.parent()
    weights = PR.term_order().weights()
    r = weights.count(1)
    m = max(weights) + 1
    q = Q.degree()

    if low_memory:
        return _step_upper_bound_low_mem(r, m, q, generating_function)
    return _step_upper_bound_internal(r, m, q, generating_function)


def first_solvable_system(Q, min_step=None):
    r"""
    Return the matrix in the lowest step such that there exists a nontrivial
    abnormal covector admitting the factor Q.

    INPUT:

    - ``Q`` -- a polynomial to search for as an abnormal factor
    - ``min_step`` -- the lowest nilpotency step to check for solutions

    OUTPUT:

    A pair (params, A), where

    - ``params`` -- a tuple (L,m,s) of data related to the matrix A
    - ``A`` -- the matrix of the first solvable system

    The parameters are

    - ``L`` -- the free Lie algebra quotient
    - ``m`` -- the lowest degree of covectors in the system
    - ``s`` -- the nilpotency step where a solution is found
    """
    PR = Q.parent()
    weights = PR.term_order().weights()
    r = weights.count(1)
    m = max(weights) + 1
    s = m + Q.degree()
    if min_step is not None:
        s = max(min_step, s)

    R = PR.base_ring()
    L = NotAFreeLieAlgebra(R, ['X_%d' % (k + 1) for k in range(r)])
    L = L.HallQuotient(m)

    systems = {}
    while True:
        verbose("searching for solutions in step %d" % s)
        A = abnormal_factor_system(L, Q, r, m, s, homogeneous=False)
        systems[s] = A

        if A.rank() < len(A.columns()):
            verbose("first solution exists in step %d" % s)
            break

        s += 1

    return (L, m, s), A


def extract_a_solution(L, m, s, A, homogeneous=False, sol_type="any"):
    r"""
    Return a single solution from the kernel of the matrix A.

    INPUT:

    - ``L`` -- a free Lie algebra quotient
    - ``m`` -- the lowest degree of covectors in the matrix
    - ``s`` -- the highest degree of covectors in the matrix
    - ``A`` -- a matrix of the form given by :func:`abnormal_factor_system`
    - ``homogeneous`` -- a boolean; whether the solved polynomial is homogeneous
    - ``sol_type`` -- one of "any", "all", "nice"

    With sol_type="any", the output is the solution from the first nonpivot
    column of the echelon form of A.
    With sol_type="all", the output is the basis of the kernel given by
    the nonpivot columns of the echelon form of A.
    With sol_type="nice", the output is a linear combination that attempts
    to decrease the number of nonzero coefficients and their magnitudes.

    OUTPUT:

    A dictionary {X: c} describing a covector. Each `X` is an element of the
    Lie algebra `L` and each `c` is the coefficient of the covector
    for the dual element of `X`.
    """
    # extract names of the covector variables
    if homogeneous:
        vecs = list(L.graded_basis(s))
    else:
        vecs = [X for k in range(m, s + 1) for X in L.graded_basis(k)]

    # check that the dimensions are correct
    assert A.subdivisions()[1][0] == len(vecs)
    covecdim = len(vecs)

    # compute the echelon form of A
    EA = A.echelon_form()

    def covec_to_dict(v):
        covec = {}
        for k, vk in enumerate(v):
            if vk:
                covec[vecs[k]] = vk
        return covec

    nonpivot = covecdim
    covecs = []
    while True:
        # find next non-pivot column

        while nonpivot in A.pivots():
            nonpivot += 1
        if nonpivot >= A.dimensions()[1]:
            if sol_type != "any":
                break

            # no solution found, something is wrong
            raise ValueError("something went wrong: no nonpivot column in the multiplier variables")

        if sol_type == "any":
            verbose("computing a solution")
        else:
            verbose("computing solution #%d" % (len(covecs) + 1))

        # the non-pivot columns of the echelon form define the valid solutions
        # rescale the column so that the covector has simplified coefficients
        cv = EA.column(nonpivot)[:covecdim]
        cv = cv / gcd(cv)

        if sol_type == "any":
            return covec_to_dict(cv)

        covecs.append(cv)
        nonpivot += 1

    if sol_type == "all":
        return [covec_to_dict(v) for v in covecs]

    if A.base_ring() == QQ:
        # over rationals, start with the covector with smallest coefficients
        covecs = sorted(covecs, key=lambda v: min(abs(vk) for vk in v))
    covec = covecs[0]
    covecs = covecs[1:]

    while covecs:
        other = covecs[0]
        covecs = covecs[1:]
        # eliminate the coefficient with the largest gcd
        gcdvec = [gcd(a, b) for a, b in zip(covec, other)]
        i = gcdvec.index(max(gcdvec))
        covec = covec[i] * other - other[i] * covec
        covec = covec / gcd(covec)

    # make leading coefficient positive
    i = 0
    while not covec[i]:
        i += 1
    if covec[i] < 0:
        covec = -covec

    return covec_to_dict(covec)


def ode_search(ode, odevars, sol_type="any",
               min_step=None, step_bound=True):
    r"""
    Return an abnormal covector for trajectories of the ODE through zero.

    INPUT:

    - ``ode`` -- a list of polynomials
    - ``odevars`` -- the variables of the ode
    - ``sol_type`` -- one of "any", "all", "nice"
    - ``min_step`` -- the lowest nilpotency step to check for solutions
    - ``step_bound`` -- a boolean; whether to compute an a priori bound for
      the latest nilpotency step where a solution will be found

    With sol_type="any", the output is any solution.
    With sol_type="all", the output is a basis of solutions in the lowest
    step quotient Lie algebra where a solution exists.
    With sol_type="nice", the output is a linear combination that attempts
    to decrease the number of nonzero coefficients and their magnitudes.

    OUTPUT:

    A pair (covec, aux_data), where covec is a dictionary or list of
    dictionaries describing the covector, and aux_data is a :class:`AuxData`
    containing the various objects constructed during the search.
    """
    stime = time()

    # steps 1-2: pick an orthogonal horizontal gradient and compute commutators
    F, pde = first_integral_pde(ode, odevars)
    # step 3: solve the pde
    Q = integrate_pde(F, pde)
    # step 4: compute a naive bound for when the algorithm stops
    if step_bound:
        s_max, gf = step_upper_bound_specific(Q, generating_function=True)
        verbose("upper bound for required nilpotency step is %d" % s_max)
    else:
        s_max = None
        gf = None
    # steps 5-6: increasing step one by one, check for solutions
    (L, m, s), A = first_solvable_system(Q, min_step=min_step)
    # extract one solution from the first solvable system
    covec = extract_a_solution(L, m, s, A, False, sol_type)

    etime = time()

    # write a human readable time string
    runtime = etime - stime
    if runtime >= 60:
        runtime = int(runtime)
        mins = runtime // 60
        secs = runtime % 60
        if mins >= 60:
            hours = mins // 60
            mins = mins % 60
            if hours >= 24:
                days = hours // 24
                hours = hours % 24
                timestr = "%d days %d hours %d minutes %s seconds" % (days, hours, mins, secs)
            else:
                timestr = "%d hours %d minutes %s seconds" % (hours, mins, secs)
        else:
            timestr = "%d minutes %s seconds" % (mins, secs)
    else:
        digits = ceil(-log(runtime, 10)) + 1
        if digits > 0:
            timestr = ("{:.%df} seconds" % (digits)).format(float(runtime))
        else:
            timestr = "%d seconds" % (int(runtime))
    verbose("abnormal covector found in %s" % timestr)

    return covec, AuxData(L, m, s, Q, s_max, gf)
