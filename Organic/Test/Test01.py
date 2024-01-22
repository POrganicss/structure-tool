import __init__
from Tool.Compute import Compute
from Tool.Draw import Draw

import numpy as np

y=np.arange(100)

dy=Compute.getderivative(y)

sy=Compute.getantiderivative(y)

Draw.draw_line(dy)