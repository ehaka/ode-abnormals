from sage.algebras.lie_algebras.free_lie_algebra import LieGenerator
from sage.arith.misc import factor
from sage.misc.misc import repr_lincomb


def indexword(Xs):
    if isinstance(Xs, LieGenerator):
        return str(Xs).split("_")[-1]
    return "".join([a.split("_")[-1] for a in Xs.to_word()])


def freebasis_element_to_short_word(X):
    return "X_%s" % indexword(X.leading_support())


def factorstr_withfinal(Xs, m):
    r"""
    Return a factorization string of a Lie bracket.

    INPUT:

    - ``Xs`` -- a :class:`GradedLieBracket`
    - ``m`` -- the degree of the final factor to aim for

    In the special case ``m=0``, the word is not factored.

    OUTPUT:

    A pair (f, Y), where f is the factorization string and Y is the final
    :class:`GradedLieBracket` of the factorization.
    The final factor is guaranteed to have degree at most m.
    """
    if m == 0:
        return indexword(Xs), Xs

    if getattr(Xs, '_grade', 1) <= m:
        return "(%s)" % indexword(Xs), Xs

    # extract all identical left factors
    l = Xs._left
    i = 1
    r = Xs._right
    while getattr(r, '_grade', 1) > m and r._left == l:
        i += 1
        r = r._right

    if i == 1:
        fstr = "(%s)" % (indexword(l))
    else:
        fstr = "(%s)^{%d}" % (indexword(l), i)
    tailstr, final = factorstr_withfinal(r, m)
    return fstr + tailstr, final


def covector_txt(covec, m=0, factor_coeff=False):
    r"""
    Return a txt description of the covector.

    INPUT:

    - ``covec`` -- a dictionary in the format of :meth:`extract_a_solution`
    - ``m`` -- the stopping degree of the factorization of the covector
    - ``factor_coeff`` -- (default:``False``) a boolean; whether to factor the 
      coefficients for the output

    OUTPUT:

    A text string whose rows describe the coefficients of the covector, e.g.

    (12)^{2}(2)(1)^{3}(112):-15

    means that the covector contains a summand `-15\lambda_{12122111112}`.

    The components are given as factored strings with the final term
    having degree as close to `m` as possible. The coefficients are integers.
    """

    rows = []
    for X in covec:
        Xs = X.leading_support()
        if factor_coeff:
            coeff = str(factor(covec[X]))
        else:
            coeff = str(covec[X])
        fstr, final = factorstr_withfinal(Xs, m)
        rows.append(("%s:%s" % (fstr, coeff), final, Xs))
    sortkey = lambda row: (getattr(row[2], '_grade', 1), row[1], row[2])
    sortedrows = [fstr for fstr, final, Xs in sorted(rows, key=sortkey)]

    return "\n".join(sortedrows)


def covector_latex(covec, m=0):
    r"""
    Return a latex string of the covector.

    INPUT:

    - ``covec`` -- a dictionary in the format of :meth:`extract_a_solution`
    - ``m`` -- the stopping degree of the factorization of the covector
    """
    data = []
    for X in covec:
        Xs = X.leading_support()
        fstr, final = factorstr_withfinal(Xs, m)
        data.append(("\\lambda_{%s}" % fstr, covec[X], final, Xs))
    sortkey = lambda d: (getattr(d[3], '_grade', 1), d[2], d[3])
    sortedterms = [(l, c) for l, c, f, Xs in sorted(data, key=sortkey)]

    return repr_lincomb(sortedterms, is_latex=True)

