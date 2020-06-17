from .print_utilities import freebasis_element_to_short_word

from copy import copy
from itertools import combinations
from sage.algebras.lie_algebras.lie_algebra import LieAlgebra
from sage.algebras.lie_algebras.free_lie_algebra import (GradedLieBracket,
                                                         LieGenerator)
from sage.misc.misc import verbose
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.rational_field import QQ


def first_integral_pde(ode, odevars):
    r"""
    Transform a polynomial ODE system in `\mathbb{R}^r` to a PDE for a
    first integral in a free Lie algebra with `r` generators.

    INPUT:

    - ``ode`` -- a list of polynomials giving the derivatives of each variable
    - ``odevars`` -- a list of the variables

    OUTPUT:

    A pair ``(F, pde)``, where

    - ``F`` -- a free Lie algebra in a Hall basis
    - ``pde`` -- a dictionary whose keys are supports of basis elements of ``F``
      and values are the partial derivatives with respect to the basis element
    """
    if any(eq.is_zero() for eq in ode):
        raise ValueError("ode must not contain a trivial equation")

    # define free Lie algebra with number of generators the variables
    PR = ode[0].parent()
    R = PR.base_ring()
    r = len(odevars)
    F = LieAlgebra(R, ['X_%d' % (k + 1) for k in range(r)]).Hall()

    operators = {X.leading_support(): lambda P, v = ov: P.derivative(v)
                 for X, ov in zip(F.gens(), odevars)}

    # choose two equations of minimal degree from the ode and define
    # the first integral to have horizontal derivatives orthogonal to those two
    degs = [(ode[i].degree(), i) for i in range(len(ode))]
    (ideg, i), (jdeg, j) = sorted(degs)[:2]
    maxdeg = max(ideg, jdeg)
    pde = {X.leading_support(): PR.zero() for X in F.gens()}
    pde[F.gens()[i].leading_support()] = ode[j]
    pde[F.gens()[j].leading_support()] = -ode[i]

    verbose("computing commutators from the partial derivatives:")
    for xk in pde:
        verbose("    %s Q = %s" % (xk, pde[xk]))
    # compute all commutators up to maxdeg+1 or until all derivatives are zero
    zeromap = lambda P: PR.zero()
    for k in range(2, 2 + maxdeg):
        allzero = True

        for X in F.graded_basis(k):
            # compute the commutator derivative [a,b]Q = a(bQ)-b(aQ)
            ab = X.leading_support()
            a = ab._left
            b = ab._right
            opa = operators.get(a, zeromap)
            opb = operators.get(b, zeromap)
            pde_ab = opa(pde[b]) - opb(pde[a])
            if pde_ab:
                verbose("    %s Q = %s" % (ab, pde_ab))
            pde[ab] = pde_ab
            allzero = allzero and pde_ab.is_zero()

        if allzero:
            # remove unnecessary zeros from pde
            for X in F.graded_basis(k):
                ab = X.leading_support()
                del pde[ab]
            break

    return F, pde

def integrate_pde(F, pde):
    """
    Integrate a pde for a first integral of an ode in the free Lie algebra.

    INPUT

    - ``F`` -- free Lie algebra whose element supports are the keys of ``pde``
    - ``pde`` -- the pde dictionary output by :func:`first_integral_pde`
    """
    # compute rank and step for the ambient free nilpotent Lie algebra
    rank = len([a for a in pde if isinstance(a, LieGenerator)])
    if len(pde) == rank:
        # everything is horizontal
        step = 1
    else:
        step = max(a._grade for a in pde if isinstance(a, GradedLieBracket))
    pdevars = [a for a in pde if pde[a]]
    PDE_PR = list(pde.values())[0].parent()
    R = PDE_PR.base_ring()

    # define a free nilpotent Lie algebra of rank r step s in the Hall basis
    # by extracting structure coefficients from the free Lie algebra
    verbose("extracting structure coefficients for a free nilpotent Lie algebra or rank %d and step %d" % (rank, step))
    shortnames = {}
    for X in F.basis():
        s = X.leading_support()
        if isinstance(s, GradedLieBracket) and s._grade > step:
            break
        shortnames[s] = freebasis_element_to_short_word(X)
    sc = {}
    for s in range(2, step + 1):
        for k in range(1, (s + 1) // 2):
            # brackets [V_k, V_{s-k}]
            for X in F.graded_basis(k):
                Xs = X.leading_support()
                for Y in F.graded_basis(s - k):
                    Ys = Y.leading_support()
                    Z = X.bracket(Y)
                    if Z:
                        Zdict = {shortnames[W]:c for W, c
                                 in Z.monomial_coefficients().items()}
                        sc[(shortnames[Xs], shortnames[Ys])] = Zdict
        # if s is even, consider also brackets [V_{s/2},V_{s/2}]
        if s % 2 == 0:
            for X, Y in combinations(F.graded_basis(s // 2), 2):
                Xs = X.leading_support()
                Ys = Y.leading_support()
                Z = X.bracket(Y)
                if Z:
                    Zdict = {shortnames[W]:c for W, c
                             in Z.monomial_coefficients().items()}
                    sc[(shortnames[Xs], shortnames[Ys])] = Zdict
    L = LieAlgebra(QQ, sc, names=list(shortnames.values()), nilpotent=True)

    # define a polynomial ring in the variables of the Lie algebra
    weights = [1] * rank + [s._grade for s in shortnames
                          if isinstance(s, GradedLieBracket)]
    PR = PolynomialRing(R, [s.lower() for s in shortnames.values()],
                        order=TermOrder('wdegrevlex', weights))

    # compute derivations on the polynomial ring from
    # exp-2 coords for the free nilpotent Lie algebra
    verbose("computing left invariant vector fields in exp2 coordinates")
    G = L.lie_group()
    exp2 = G.chart_exp2()
    G.set_default_chart(exp2)
    G.set_default_frame(exp2.frame())
    dLx = G._dLx()
    operators = {}
    vfsubs = dict(zip(exp2[:], PR.gens()))

    def derivation(P, vf):
        XP = PR.zero()
        for vk, xk in zip(vf, PR.gens()):
            XP += PR(vk.subs(vfsubs)) * P.derivative(xk)
        return XP

    for vf, xk in zip(dLx.columns(), PR.gens()):
        operators[xk] = lambda P, vf = vf: derivation(P, vf)

    # integrate the pde in each variable in reverse Hall order
    Q = PR.zero()
    varsubs = {ak:bk for ak, bk in zip(PDE_PR.gens(), PR.gens())}
    pde = {xk:PR(pde[Xk].subs(varsubs)) for xk, Xk in zip(PR.gens(), pde)}
    remaining = copy(pde)

    verbose("integrating the PDE")
    for xk in reversed(PR.gens()):
        # integrate in xk
        R = remaining[xk].integral(xk)
        if not R:
            continue
        Q += R
        for yk in remaining:
            remaining[yk] -= operators[yk](R)

        i = PR.gens().index(xk)
        remainingvarstr = ""
        if i > 0:
            remainingvarstr += " + R("
            remainingvarstr += ",".join(str(yk) for yk in PR.gens()[:i])
            remainingvarstr += ")"
        verbose("  Q = %s%s" % (Q, remainingvarstr))

    # verify correctness of solution
    verbose("verifying the solution")
    for xk in pde:
        assert operators[xk](Q) == pde[xk]

    return Q
