from .print_utilities import freebasis_element_to_short_word

from copy import deepcopy
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
from sage.functions.other import factorial
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc import verbose
from sage.rings.rational_field import QQ


@cached_method
def abnormal_polynomials(L, r, m, s):
    r"""
    Return the highest degree coefficients for abnormal polynomials in step `s`
    of layer `m` in the Lie algebra `L`.

    INPUT:

    - ``L`` -- the Lie algebra to do the computations in
    - ``r`` -- the rank of the free Lie algebra
    - ``m`` -- the layer whose abnormal polynomials are computed
    - ``s`` -- the step of coefficients to compute

    OUTPUT:

    A dictionary {X: {mon: coeff}, where

    - ``X`` -- a Hall basis element of layer `m`
    - ``mon`` -- a tuple `(i_1,...,i_n)` describing the monomial
      `x_1^{i_1}...x_n^{i_n}`
    - ``coeff`` -- an element of a :class:`HallQuotient` Lie algebra whose dual
      element would be the coefficient of the monomial in the polynomial `P_X`
    """
    elems = sum((list(L.graded_basis(k)) for k in range(1, m)), [])
    weights = sum(([k] * len(L.graded_basis(k)) for k in range(1, m)), [])
    P = {}
    mons = list(reversed(list(WeightedIntegerVectors(s - m, weights))))
    verbose("computing abnormal polynomials:")
    for X in L.graded_basis(m):
        PX = {}
        pname = freebasis_element_to_short_word(X).replace("X", "P")
        verbose("  %s deg %d coefficients" % (pname, s - m))
        # compute coefficients for the abnormal polynomial P_X
        for mon in mons:
            # compute coefficient of the monomial mon
            adx = X
            mult = QQ(1)
            for i in mon:
                mult = mult / factorial(i)
            for rep, Xi in zip(mon, elems):
                for k in range(rep):
                    adx = Xi.bracket(adx)
            PX[tuple(mon)] = mult * adx
        P[X] = PX

    return P


def abnormal_polynomials_with_covec(PR, covec, m):
    r"""
    Return abnormal polynomials with a given covector.

    INPUT:

    - ``PR`` -- the polynomial ring where the polynomial lives
    - ``covec`` -- a dictionary in the format of :func:`extract_a_solution`
    - ``m`` -- the layer of abnormal polynomials to return
    """
    # list all relevant monomials up to degree s-m-deg(Q)
    weights = PR.term_order().weights()

    # extract the quotient Lie algebra and the rank and step
    L = next(iter(covec.keys())).parent()
    r = len(L.gens())
    s = max(X.leading_support()._grade for X in covec)

    # define a function to substitute the covector
    scovec = {X.leading_support(): c for X, c in covec.items()}

    zero = L.base_ring().zero()

    def covecsubs(X):
        sub = zero
        for m, c in X.monomial_coefficients().items():
            sub += c * scovec.get(m, zero)
        return sub

    # combine the homogeneous components up to degree s with substitutions
    P = {X: zero for X in L.graded_basis(m)}
    for k in range(m, s + 1):
        Pk = abnormal_polynomials(L, r, m, k)
        for X in P:
            for mon in Pk[X]:
                P[X] += covecsubs(Pk[X][mon]) * PR.monomial(*mon)
    return P


def abnormal_factor_system(L, Q, r, m, s, homogeneous=True):
    r"""
    Return a matrix for the linear system to find `Q` as a common factor
    in the abnormal polynomials of layer `m` in rank `r`.

    INPUT:

    - ``L`` -- the Lie algebra to do the computations in
    - ``Q`` -- the polynomial to look for as a common factor
    - ``r`` -- the rank of the Lie algebra
    - ``m`` -- the layer of abnormal polynomials to search
    - ``s`` -- the step to search in
    - ``homogeneous`` -- a boolean; if ``True`` only the highest order component
      of the matrix is computed. Otherwise the full system is computed.

    OUTPUT:

    A pair (L,A), where `L` is the ambient free Lie algebra quotient and
    `A` is a block matrix of coefficients of the linear system. The variables
    are ordered so that the first block has `\dim f_r^{s}` columns describing
    the covector and the rest are the auxiliary multiplier polynomial variables.
    """
    PR = Q.parent()
    weights = PR.term_order().weights()

    # list all relevant monomials up to degree s-m-deg(Q)
    def monomial_list(deg):
        return reversed(list(WeightedIntegerVectors(deg, weights)))

    if homogeneous:
        mons = [tuple(mon) for mon in monomial_list(s - m)]
        smons = [PR.monomial(*mon) for mon in monomial_list(s - m - Q.degree())]
        # extract the maximal degree monomials of the polynomial Q
        qmons = [qm for qm in Q.monomials() if qm.degree() == Q.degree()]
        P = abnormal_polynomials(L, r, m, s)
        supp_basis = [X.leading_support() for X in L.graded_basis(s)]
    else:
        mons = [tuple(mon) for k in range(s - m + 1)
                           for mon in monomial_list(k)]
        smons = [PR.monomial(*mon) for k in range(s - m - Q.degree() + 1)
                                   for mon in monomial_list(k)]
        qmons = Q.monomials()
        P = None
        for k in range(m, s + 1):
            Pk = abnormal_polynomials(L, r, m, k)
            if not P:
                # make a copy to not ruin cache
                P = deepcopy(Pk)
                continue
            for X, PX in Pk.items():
                for Y, c in PX.items():
                    P[X][Y] = c
        supp_basis = [X.leading_support() for k in range(m, s + 1)
                      for X in L.graded_basis(k)]

    def lie_elem_to_dict(X):
        mc = X.monomial_coefficients()
        return {supp_basis.index(m): c for m, c in mc.items()}

    # form the linear system as a sparse matrix
    typestr = "homogeneous" if homogeneous else "full"
    verbose("forming %s linear system in step %d" % (typestr, s))
    covec_block = {}
    aux_blocks = {X: {} for X in P}
    i = 0
    for mon in mons:
        pmon = PR.monomial(*mon)

        # S parameters are all coefficients of monomials
        sc = {}
        for qmon in qmons:
            q, r = pmon.quo_rem(qmon)
            if r.is_zero() and q.degree() <= s - m - Q.degree():
                c = Q.monomial_coefficient(qmon)
                if q in smons:
                    sc[smons.index(q)] = -c

        for X, PX in P.items():
            # compute the coefficient of mon in P_X - S_X * Q
            # as a row vector of the coefficients of the P and S params

            # P parameters are the covector
            pc = PX.get(tuple(mon), L.zero())
            for j, c in lie_elem_to_dict(pc).items():
                covec_block[(i, j)] = c

            for j, c in sc.items():
                aux_blocks[X][(i, j)] = c

            i += 1

    numrows = i
    verbose("%d equations, %d covector vars, %d x %d auxiliary vars" % (numrows,
            len(supp_basis), len(P), len(smons)))
    R = PR.base_ring()
    covecmat = matrix(R, numrows, len(supp_basis), covec_block)
    auxmats = [matrix(R, numrows, len(smons), d) for d in aux_blocks.values()]
    return matrix.block(R, [covecmat] + auxmats, ncols=1 + len(auxmats))
