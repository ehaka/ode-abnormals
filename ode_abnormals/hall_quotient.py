from six import iteritems
from sage.algebras.lie_algebras.free_lie_algebra import (FreeLieAlgebra,
                                                         FreeLieBasis_abstract,
                                                         GradedLieBracket,
                                                         LieGenerator)
from sage.misc.cachefunc import cached_method
from sage.structure.indexed_generators import standardize_names_index_set


class NotAFreeLieAlgebra(FreeLieAlgebra):

    @staticmethod
    def __classcall_private__(cls, R, names, index_set=None):
        """
        Normalize input to ensure a unique representation.
        """
        names, index_set = standardize_names_index_set(names, index_set)
        return super(NotAFreeLieAlgebra, cls).__classcall__(cls, R, names,
                                                            index_set)

    class HallQuotient(FreeLieAlgebra.Hall):
        """
        Helper class for Hall trees in the quotient f/[f_k, f_k],
        where f_k is the k:th term of the lower central series.

        .. WARNING::

            By an abuse of existing code, the construction is treated as a basis
            for a free Lie algebra, which may cause problems with careless use.
        """

        def __init__(self, lie, q):
            FreeLieBasis_abstract.__init__(self, lie,
                                           "Hall quotient of degree %d" % q)
            self._maxdeg = q

        # override without implementation to avoid incorrect data
        def graded_dimension(self, k):
            return NotImplemented

        @cached_method
        def _rewrite_bracket(self, l, r):
            if isinstance(l, GradedLieBracket) and isinstance(r, GradedLieBracket):
                if l._grade >= self._maxdeg and r._grade >= self._maxdeg:
                    return {}

            ####################################################################
            # The rest of the method is copied verbatim from
            # sage.algebras.lie_algebras.free_lie_algebra.FreeLieAlgebra.Hall
            # Calling super()._rewrite_bracket instead of copying
            # the code here breaks the method caching
            ####################################################################

            if not isinstance(r, GradedLieBracket) or r._left <= l:
                # Compute the grade of the new element
                grade = 0
                # If not a graded Lie bracket, it must be a generator so the grade is 1
                if isinstance(l, GradedLieBracket):
                    grade += l._grade
                else:
                    grade += 1
                if isinstance(r, GradedLieBracket):
                    grade += r._grade
                else:
                    grade += 1
                return {GradedLieBracket(l, r, grade): self.base_ring().one()}

            ret = {}

            # Rewrite [a, [b, c]] = [b, [a, c]] + [[a, b], c] with a < b < c
            # Compute the left summand
            for m, inner_coeff in iteritems(self._rewrite_bracket(l, r._right)):
                if r._left == m:
                    continue
                elif r._left < m:
                    x, y = r._left, m
                else:  # r._left > m
                    x, y = m, r._left
                    inner_coeff = -inner_coeff
                for b_elt, coeff in iteritems(self._rewrite_bracket(x, y)):
                    ret[b_elt] = ret.get(b_elt, 0) + coeff * inner_coeff

            # Compute the right summand
            for m, inner_coeff in iteritems(self._rewrite_bracket(l, r._left)):
                if m == r._right:
                    continue
                elif m < r._right:
                    x, y = m, r._right
                else:  # m > r._right
                    x, y = r._right, m
                    inner_coeff = -inner_coeff
                for b_elt, coeff in iteritems(self._rewrite_bracket(x, y)):
                    ret[b_elt] = ret.get(b_elt, 0) + coeff * inner_coeff

            return ret

        @cached_method
        def _generate_hall_set(self, k):
            """
            Generate the Hall set of grade ``k``.
            """
            m = self._maxdeg
            if k >= 2 * m:
                # adapt the construction of Hall._generate_hall_set to skip [a,b]
                # when both a and b have grade>=self._maxdegree
                ret = [GradedLieBracket(a, b, k) for i in range(1, m)
                       for a in self._generate_hall_set(i)
                       for b in self._generate_hall_set(k - i)
                       if b._left <= a]

                return tuple(sorted(ret))

            ####################################################################
            # The rest of the method is copied verbatim from
            # sage.algebras.lie_algebras.free_lie_algebra.FreeLieAlgebra.Hall
            # Calling super()._generate_hall_set instead of copying
            # the code here breaks the method caching
            ####################################################################

            if k <= 0:
                return ()
            if k == 1:
                return tuple(map(LieGenerator, self.variable_names(),
                                 range(len(self.variable_names()))))
            if k == 2:
                basis = self._generate_hall_set(1)
                ret = []
                for i, a in enumerate(basis):
                    for b in basis[i + 1:]:
                        ret.append(GradedLieBracket(a, b, 2))
                return tuple(ret)

            # We don't want to do the middle when we're even, so we add 1 and
            #   take the floor after dividing by 2.
            ret = [GradedLieBracket(a, b, k) for i in range(1, (k + 1) // 2)
                   for a in self._generate_hall_set(i)
                   for b in self._generate_hall_set(k - i)
                   if b._left <= a]

            # Special case for when k = 4, we get the pairs [[a, b], [x, y]]
            #    where a,b,x,y are all grade 1 elements. Thus if we take
            #    [a, b] < [x, y], we will always be in the Hall set.
            if k == 4:
                basis = self._generate_hall_set(2)
                for i, a in enumerate(basis):
                    for b in basis[i + 1:]:
                        ret.append(GradedLieBracket(a, b, k))
            # Do the middle case when we are even and k > 4
            elif k % 2 == 0:
                basis = self._generate_hall_set(k // 2)  # grade >= 2
                for i, a in enumerate(basis):
                    for b in basis[i + 1:]:
                        if b._left <= a:
                            ret.append(GradedLieBracket(a, b, k))

            # We sort the returned tuple in order to make computing the higher
            #   graded parts of the Hall set easier.
            return tuple(sorted(ret))
