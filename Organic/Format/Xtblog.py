import os
import sys

from Log import Log
path = os.getcwd()
path = path.replace("\\", "/")+'/Organic'
sys.path.append(path+'/Format')
sys.path.append(path+'/Tool')
sys.path.append(path+'/Applications')
sys.path.append(path+'/Functions')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.widgets import Slider
import mplcursors
from Compute import *
from Datatransmission import *
import numpy as np
from Xyz import *
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

    # 获取log中初始参数
    def getparameters(log):
        rline = log.splitlines()
        paramenters = {}
        for j in range(len(rline)):
            if '%mem' in rline[j]:
                paramenters['mem'] = rline[j].split("=")[1].strip()
            elif '%cpu' in rline[j]:
                paramenters['cpu'] = rline[j].split("=")[1].strip()
            elif '#' in rline[j]:
                paramenters['code'] = rline[j].strip()
            elif 'Multiplicity =' in rline[j] and ' Charge = ' in rline[j]:
                paramenters['charge'] = rline[j].split()[2].strip()
                paramenters['spin'] = rline[j].split()[5].strip()
            elif ' Input=' in rline[j]:
                paramenters['name'] = rline[j].split("=")[1].strip()[:-4]
        return paramenters

    # 获取gaussian运行时的参数
    def getopearationinformations(log):
        Information = {}
        Information['SPE'] = []  # 单点能

        Information['MaxForce'] = []  # 最大受力
        Information['RMSForce'] = []  # 最大受力的根均方

        Information['MaxDisplacement'] = []  # 最大位移
        Information['RMSDisplacement'] = []  # 最大位移的根均方

        Information['Info'] = []  # 最大位移的根均方

        log_lines = log.splitlines()

        for line in log_lines:
            if 'SCF Done' in line:
                Information['SPE'].append(float(line.split()[4]))
            if 'Maximum Force' in line:
                Information['MaxForce'].append(float(line.split()[2]))
            if 'RMS     Force' in line:
                Information['RMSForce'].append(float(line.split()[2]))
            if 'Maximum Displacement' in line:
                Information['MaxDisplacement'].append(float(line.split()[2]))
            if 'RMS     Displacement' in line:
                Information['RMSDisplacement'].append(float(line.split()[2]))
            if 'Internal  Forces' in line:
                Information['Info'].append(float(line.split()[3]))
        return Information

    # 获取gaussian运行时的某参数的变化图像
    def getdraw(log, babel):
        # 从日志获取操作信息
        Info = Log.getopearationinformations(log)

        # 使用提供的键提取特定数据
        data = Info[babel]
        if len(data) == 0:
            print('没有该参数')
        else:
            Compute.draw(data)
  
data = File.getdata('C:\\Users\\10282\\OneDrive\\桌面\\xtbscan.log')
xyzs,en=Xtblog.toxyz(data)
x = Xyz.getparameter(xyzs, 8, 5, 1, 4)
Compute.draws(x, en)

