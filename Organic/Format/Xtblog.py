import __init__
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
import mplcursors

from Format.Log import Log
from Format.Xyz import *
from Tool.Compute import *
from Tool.Datatransmission import *
import numpy as np

class Xtblog:
    def toxyz(log):
        xyzs = []
        current_xyz = []
        energy=[]
        for line in log.splitlines():
            if 'energy:' in line:
                energy.append(float(line.split()[1]))
                continue
            
            string = line.split()
            if len(string) == 4 and not Verify.isnums(string[0]) and Verify.isnums(string[1:4]):
                current_xyz.append(line + '\n')
                
            elif len(current_xyz) != 0:
                xyz = Xyz.initialize(' '.join(current_xyz))
                xyzs.append(xyz)
                current_xyz = []
        if current_xyz != []:
            xyz = Xyz.initialize(' '.join(current_xyz))
            xyzs.append(xyz)
        return xyzs,energy

  
data = File.getdata('C:\\Users\\10282\\OneDrive\\桌面\\xtbscan.log')
xyzs,en=Xtblog.toxyz(data)
x = Xyz.getparameter(xyzs, 8, 5, 1, 4)
Compute.draws(x, en)

