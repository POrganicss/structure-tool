import __init__
from tkinter import messagebox
from Tool.Code import Code
from Tool.Verify import Verify

class Gjf:
    #所有的参数：name,mem,nproc,code,mixset
    #其中基组的存储形式sets=[[['C'，'H','O'],'6-31g'],[['Pd'],'SDD']]
    #读取形式：elements=[['C'，'H','O'],['Pd']],mixset=['6-31g','SDD']
    
    # 从gjf中获取xyz坐标文件
    @staticmethod
    def toxyz(gjf):
        lines = gjf.splitlines()  # 将输入数据拆分为行
        xyzs = []
        xyz = []  
        for line in lines:
            string = line.split()
            if len(string) == 4 and not Verify.isnum(string[0]) and Verify.isnums(string[1:4]):
                xyz.append(line)
            elif len(xyz) != 0:
                xyzs.append('\n'.join(xyz))
                xyz = []
        if len(xyzs) <= 1:
            return xyzs[0]
        else:
            return xyzs

    # 从gjf中获取除坐标外的所有参数
    def getcode(gjf):
        
        def get_functional(chunk):
            readline=chunk.splitlines()
            paramenters = {}
            for line in readline:
                line = line.strip()
                if '%chk' in line:
                    paramenters['name'] = line.split("=")[1].split(".")[0]
                elif '%mem' in line:
                    paramenters['mem'] = line.split("=")[1]
                elif '%nproc' in line:
                    paramenters['nproc'] = line.split("=")[1]
                elif '%cpu' in line:
                    paramenters['cpu'] = line.split("=")[1]
                elif '#' in line:
                    paramenters['code'] = line
                elif len(line.split()) == 2 and Verify.isnums(line.split()):
                    paramenters['charge'] = line.split()[0]
                    paramenters['spin'] = line.split()[1]
            return paramenters
     
        def get_mixset(chunk):
            lines=chunk.strip().splitlines()
            
            elements=[]#存储类似elements=[['C'，'H','O'],['Pd']]
            MIXSET=[]#存储类似mixset=['6-31g','SDD']
            
            elements.append(lines[0].split())
            elements[0].remove('0')
            MIXSET.append(lines[1])
            
            for i in range(1, int(len(lines)/3)):
                element = lines[i*3].split()
                element.remove('0')
                elements.append(element)
                MIXSET.append(lines[i*3+1])
            
            return {'mixset':[[elements[0],MIXSET[0]],[elements[1],MIXSET[1]]]}
        
        chunk=gjf.replace('\n\n\n','\n\n').split('\n\n')
        Header=chunk[0]+'\n'+chunk[2].splitlines()[0]+'\n'
        Footer=chunk[3]
        para = {}
        para.update(get_functional(Header))
        para.update(get_mixset(Footer))
        return para
