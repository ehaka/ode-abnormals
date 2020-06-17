# Import the code from the src folder without installation
# by temporarily adding the src folder to path.
import sys
import pathlib
path = pathlib.Path().absolute().parent
sys.path.append(str(path))

from ode_abnormals.ode_search import step_upper_bound
from time import time

# Compute the a priori bound for abnormality step with two methods.
r = 4
d = 4
at = time()
s = step_upper_bound(r, d, low_memory=True)
bt = time()
print("s(%d,%d) = %d, computed with low mem in %.1f seconds\n" % (r, d, s,
                                                                  bt - at))

# Cheat series precision to see what the most efficient computation could be.
set_series_precision(s + 1)
at = time()
s = step_upper_bound(r, d, low_memory=False)
bt = time()
print("s(%d,%d) = %d, computed in %.1f seconds\n" % (r, d, s, bt - at))
