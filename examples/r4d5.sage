from time import time
from collections import deque
from datetime import datetime

curtime = datetime.now()
print("Computation started at %s"%(curtime.strftime('%H:%M:%S')))

# Compute the a priori bound for abnormality step
r = 4
d = 5
m = d + 2
q = d + 1

class LoopData:
    def __init__(self, s, coeffs, cumul, multipliers):
        self.s = s
        self.coeffs = coeffs
        self.multipliers = multipliers
        self.cumul = cumul

savefile = 'loopdata.dump'
try:
    with open(savefile,'rb') as f:
        ld = loads(f.read())
except IOError:
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
            outstr = "s(%d,%d) = %d, computed with low mem in %.1f seconds\n" % (r, d, s, bt - at)
            with open('r4d5.txt','w') as f:
                f.write(outstr)
            print(outstr)

    # the rest of the coefficients are defined by a recurrence relation
    multipliers = [-b.monomial_coefficient(t ** (bd - k)) for k in range(bd)]

    ld = LoopData(s,coeffs, cumul ,multipliers)
    data = dumps(ld)
    with open(savefile,'wb') as f:
        f.write(data)


from time import time
at = time()

saveinterval = 499999
coeffs = ld.coeffs
s = ld.s
multipliers = ld.multipliers
cumul = ld.cumul
prog = s %(saveinterval+1)

print("Starting from %d"%s)


while cumul <= 0:
    c_next = sum(c * m for c, m in zip(coeffs, multipliers))
    cumul += c_next
    s += 1
    prog+=1
    coeffs.append(c_next)
    coeffs.popleft()
    
    if prog>saveinterval:
        curtime = datetime.now()
        timestr = curtime.strftime('%H:%M:%S')
        print("Step %d, recurrence runtime %d seconds at %s"%(s,int(time()-at),timestr))
        ld = LoopData(s,coeffs, cumul,multipliers)
        data = dumps(ld)
        with open(savefile,'wb') as f:
            f.write(data)
        prog=0

outstr = "s(%d,%d) = %d, computed with low mem in %.1f seconds\n" % (r, d, s, time() - at)
print(outstr)
