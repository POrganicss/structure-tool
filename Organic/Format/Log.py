import __init__

from Xyz import *
import matplotlib.pyplot as plt
import numpy as np
from Tool.Datatransmission import *
from Tool.Compute import *
from Tool.Draw import Draw
class Log:
    elements = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',
            7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
            13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
            18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
                23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co',
            28: 'Ni', 29: 'Cu', 30: 'Zn', 33: 'As', 35: 'Br',
            44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 49: 'In',
            53: 'I', 55: 'Cs', 74: 'W', 75: 'Re', 77: 'Ir',
            78: 'Pt', 79: 'Au', 82: 'Pb'}
    #获取log中所有标准坐标文件
    
    def getxyzs(log):
        xyzs=[]
        current_xyz = []
        is_reading_coordinates = False
        delay_counter=0
        for line in log.splitlines():
            if 'Standard orientation:' in line:
                is_reading_coordinates = True
                continue
            if is_reading_coordinates:
                if delay_counter < 4:
                    delay_counter += 1
                else:
                    string = line.split()
                    if len(string) == 6 and Verify.isnums(string):
                        element = Log.elements[int(line.split()[1])]
                        current_xyz.append(element + line[30:-1] + '\n')
                    elif len(current_xyz) != 0:
                        xyz = Xyz.initialize(' '.join(current_xyz))
                        xyzs.append(xyz)
                        current_xyz = []
                        is_reading_coordinates = False
                        delay_counter = 0
        del xyzs[0]
        return xyzs
    
    def toxyz(log,index=-1):
        return Log.getxyzs(log)[index]
    
    def getcode(log):
        lines = log.splitlines()
        para = {}
        for line in lines:
            if '%mem' in line:
                para['mem']=line.split("=")[1].strip()
            elif '%cpu' in line:
                para['cpu']=line.split("=")[1].strip()
            elif '%nproc' in line or '%nprocshared' in line:
                    para['nproc']=line.split("=")[1].strip()
            elif '#' in line:
                para['code']=line.strip()
            elif 'Multiplicity =' in line and ' Charge = ' in line :
                para['charge'] = line.split()[2].strip()
                para['spin'] = line.split()[5].strip()
            elif ' Input=' in line:
                para['name']=line.split("=")[1].strip()[:-4]
             
    
                
    #获取gaussian运行时的参数
    def getopearationinformations(log):
        Information={}
        
        Information['SPE']=[]             #单点能
        
        Information['MaxForce']=[]        #最大受力
        Information['RMSForce']=[]        #最大受力的根均方
        
        Information['MaxDisplacement']=[] #最大位移
        Information['RMSDisplacement']=[] #最大位移的根均方
        
        Information['Force'] = []  # 最大位移的根均方
        
        log_lines=log.splitlines()
        
        precision_map = {
            0.01: 'loose',
            0.0018: 'default',
            0.00006: 'tight',
            0.000006: 'very tight'
        }
        
        precision_maps = {
            'loose': {
                'MaxForce': 0.0025,
                'RMSForce': 0.001667,
                'MaxDisplacement': 0.01,
                'RMSDisplacement': 0.006667,
                'Force': 0.0025
            },
            'default': {
                'MaxForce': 0.00045,
                'Force': 0.00045,
                'RMSForce': 0.0003,
                'MaxDisplacement': 0.0018,
                'RMSDisplacement': 0.0012
            },
            'tight': {
                'MaxForce': 0.000015,
                'Force': 0.000015,
                'RMSForce': 0.00001,
                'MaxDisplacement': 0.00006,
                'RMSDisplacement': 0.00004
            },
            'very tight': {
                'MaxForce': 0.000002,
                'Force': 0.000002,
                'RMSForce': 0.000001,
                'MaxDisplacement': 0.000006,
                'RMSDisplacement': 0.000004
            }
        }

        accuracy=''
        for line in log_lines:
            words = line.split()
            if 'SCF Done' in line:
                Information['SPE'].append(float(words[4]))
            elif 'Maximum Force' in line:
                Information['MaxForce'].append(float(words[2]))
            elif 'RMS     Force' in line:
                Information['RMSForce'].append(float(words[2]))
            elif 'Maximum Displacement' in line:
                Information['MaxDisplacement'].append(float(words[2]))
                if accuracy =='':
                    accuracy = precision_map[float(words[3])]
            elif 'RMS     Displacement' in line:
                Information['RMSDisplacement'].append(float(words[2]))
            elif 'Internal  Forces' in line:
                Information['Force'].append(float(words[3]))
        
        Information['limit'] = precision_maps[accuracy]
        return Information

    #获取gaussian运行时的某参数的变化图像
    
    def getdraw(log,babel):
        # 从日志获取操作信息
        Info = Log.getopearationinformations(log)
        data = Info[babel]
        if len(data)==0:
            print('没有该参数')
        else:
            limit_value = Info['limit'][babel]
            Draw.draw_line('', Info[babel], limit=[limit_value, 10])
    
    def getall(file):
        log = File.getdata(file)
        Information = Log.getopearationinformations(log)
        xyzs = Log.getxyzs(log)
        return Information,xyzs
# SPE：单点能
# MaxForce：最大受力；Force：最大力；RMSForce：最大受力的根均方
# MaxDisplacement：最大位移；RMSDisplacement：最大位移的根均方
# limit_value = Info['limit'][babel]


# Information, xyzs = Log.getall(
#     'C:\\Users\\10282\\OneDrive\\桌面\\Benzene02_Scan.log')
# B3_13=Xyz.getparameter(xyzs,2,13)
# Di12_3_13_15 = Xyz.getparameter(xyzs, 12, 3,13,15)

# # Draw.draw_line(Information['MaxDisplacement'], limit=[
# #                Information['limit']['MaxDisplacement'], 10])
# filename='C:\\Users\\10282\\OneDrive\\桌面\\data02.xlsx'
# #File.toxexcels(filename, B5_21, Di9_5_21_23, Information['Force'])
# data=Compute.renormalization(B3_13, Di12_3_13_15, Information['Force'])
# File.tofile(filename, data)


# # Draw.draw_line(Di12_3_13_15, Information['Force'], labels={
# #                'xlabel': 'D1_3', 'ylabel': 'MaxForce', 'zlabel': 'SPE', 'label': 'Draws', 'title': 'title', 'marker': 'x'})
# # # #Draw.draw_line()

