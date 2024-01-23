import __init__
from Tool.Compute import Compute
from Tool.Draw import Draw

import numpy as np

y=np.arange(10)
print(y)
dy=Compute.getderivative(y)
print(dy)

# sy=Compute.getantiderivative(y)

# Draw.draw_line(dy)